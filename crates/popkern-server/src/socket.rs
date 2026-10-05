//! The session WebSocket (`./api/session`): the worker protocol over a
//! socket, answered by a native session.
//!
//! One thread owns the socket: it reads frames with a short timeout and,
//! between reads, sends the session's replies — so frames never
//! interleave. The session runs on its own thread (requests in order,
//! like the browser worker). A `{"type": "cancel"}` frame sets the
//! session thread's cancellation flag ([`popkern_rhai::cancel`]): the
//! running script stops at its next step with the `Cancelled` class and
//! the session keeps its document — no restart. A kernel operation
//! already in progress finishes first.
//!
//! Wire format: each reply is a text frame (the protocol JSON plus
//! `attachments: [{path, kind}]`) followed by one binary frame per
//! attachment, in order; the page puts each typed array at its path.

use std::net::TcpStream;
use std::sync::atomic::Ordering;
use std::sync::mpsc;
use std::time::Duration;

use popkern_live::protocol::{Protocol, Reply};
use serde_json::{Value, json};
use tungstenite::{Message, WebSocket};

fn send(ws: &mut WebSocket<TcpStream>, reply: Reply) -> Result<(), Box<tungstenite::Error>> {
    let mut head = reply.json;
    head["attachments"] = json!(
        reply
            .attachments
            .iter()
            .map(|a| json!({ "path": a.path, "kind": a.buffer.kind() }))
            .collect::<Vec<_>>()
    );
    ws.send(Message::Text(head.to_string())).map_err(Box::new)?;
    for a in reply.attachments {
        ws.send(Message::Binary(a.buffer.bytes()))
            .map_err(Box::new)?;
    }
    Ok(())
}

/// Serve one session socket until it closes.
pub(crate) fn session(stream: TcpStream) {
    let Ok(mut ws) = tungstenite::accept(stream) else {
        return;
    };
    if ws
        .get_mut()
        .set_read_timeout(Some(Duration::from_millis(20)))
        .is_err()
    {
        return;
    }
    let (to_session, requests) = mpsc::channel::<Value>();
    let (to_socket, replies) = mpsc::channel::<Reply>();
    let (flag_out, flag_in) = mpsc::channel();
    let worker = std::thread::spawn(move || {
        let _ = flag_out.send(popkern_rhai::cancel::cancel_handle());
        let mut protocol = Protocol::new();
        for message in requests {
            if to_socket.send(protocol.handle(&message)).is_err() {
                break;
            }
        }
    });
    let Ok(cancel) = flag_in.recv() else {
        return;
    };
    'socket: loop {
        while let Ok(reply) = replies.try_recv() {
            if send(&mut ws, reply).is_err() {
                break 'socket;
            }
        }
        match ws.read() {
            Ok(Message::Text(text)) => {
                let Ok(message) = serde_json::from_str::<Value>(&text) else {
                    continue;
                };
                if message.get("type").and_then(Value::as_str) == Some("cancel") {
                    cancel.store(true, Ordering::SeqCst);
                } else if to_session.send(message).is_err() {
                    break;
                }
            }
            Ok(Message::Close(_)) => break,
            Ok(_) => {}
            Err(tungstenite::Error::Io(e))
                if matches!(
                    e.kind(),
                    std::io::ErrorKind::WouldBlock | std::io::ErrorKind::TimedOut
                ) => {}
            Err(_) => break,
        }
    }
    // A closed socket cancels what is running and ends the session.
    cancel.store(true, Ordering::SeqCst);
    drop(to_session);
    let _ = worker.join();
}
