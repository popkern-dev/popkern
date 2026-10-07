//! `popkern serve` (live-projects-and-history.md §5.7, §6.4, PK-665):
//! popkern live self-hosted. One static page for both deployments; this
//! host adds what the page asks `./host.json` for:
//!
//! - **storage = server**: projects are folders under a root you own
//!   (`<root>/<project>/popkern.toml`, `migrations/…`), read and
//!   written through `./api/projects` with revision-checked saves (a
//!   stale write answers 409 with the stored project);
//! - **compute = native** (the default; `browser` keeps the kernel in
//!   the page): `./api/session` is a WebSocket speaking the worker
//!   protocol to a native session — one session thread per
//!   connection, real cancel (`{"type": "cancel"}` stops the running
//!   script at its next step and keeps the session's document).
//!
//! - **`/mcp`**: a Model Context Protocol endpoint (JSON-RPC over HTTP
//!   POST) whose tools list, read, write and build the same projects on
//!   the native kernel, whatever `compute` the page uses.
//!
//! Local only until the security fork (PK-667) lands: binding anything
//! but a loopback address is refused ([`ServeError::NonLoopback`]).
//! There is no authentication; anyone who can reach the port can read
//! and write the projects and run scripts.

#![forbid(unsafe_code)]

mod http;
mod mcp;
mod mcp_tools;
mod socket;
pub mod store;

use std::net::{SocketAddr, TcpListener};
use std::path::PathBuf;
use std::sync::Arc;

/// Where the kernel runs for the page.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Compute {
    /// On this machine, natively, behind `./api/session`.
    Native,
    /// In the page (wasm in a Web Worker), as the public deployment.
    Browser,
}

impl Compute {
    /// `"native"` or `"browser"`.
    pub fn as_str(self) -> &'static str {
        match self {
            Compute::Native => "native",
            Compute::Browser => "browser",
        }
    }
}

/// What `popkern serve` serves.
#[derive(Debug, Clone)]
pub struct ServeOptions {
    /// The projects folder (each subfolder is a project).
    pub root: PathBuf,
    /// The page (`web/live`, with its built `pkg/`).
    pub page: PathBuf,
    /// The address to bind (loopback only).
    pub bind: SocketAddr,
    /// Where the kernel runs.
    pub compute: Compute,
}

/// Why the server could not start.
#[derive(Debug, thiserror::Error)]
pub enum ServeError {
    /// Remote access needs authentication and isolation (PK-667).
    #[error(
        "{addr} is not a loopback address: popkern serve is local only until remote \
         self-hosting (authentication, per-session isolation; ledger PK-667) lands — bind \
         127.0.0.1 or ::1"
    )]
    NonLoopback {
        /// The refused address.
        addr: SocketAddr,
    },
    /// The page folder has no `index.html`.
    #[error("{path} is not a popkern live page folder (no index.html)")]
    NoPage {
        /// The folder.
        path: PathBuf,
    },
    /// The projects root does not exist or is not a folder.
    #[error("{path} is not a folder")]
    NoRoot {
        /// The path.
        path: PathBuf,
    },
    /// The root is one project's folder, not a folder of projects.
    #[error(
        "{path} is a project folder (it holds popkern.toml): serve the folder that contains it"
    )]
    ProjectFolder {
        /// The folder.
        path: PathBuf,
    },
    /// Binding failed.
    #[error("cannot listen on {addr}: {detail}")]
    Bind {
        /// The address.
        addr: SocketAddr,
        /// The OS error.
        detail: String,
    },
}

/// Shared by every connection.
pub(crate) struct Context {
    pub(crate) options: ServeOptions,
}

/// Check `options` and bind; the listener and the context to serve.
///
/// # Errors
///
/// [`ServeError`] for a non-loopback address, a missing page or root,
/// or a failed bind.
pub fn bind(options: ServeOptions) -> Result<(TcpListener, ServeHandle), ServeError> {
    if !options.bind.ip().is_loopback() {
        return Err(ServeError::NonLoopback { addr: options.bind });
    }
    if !options.page.join("index.html").is_file() {
        return Err(ServeError::NoPage {
            path: options.page.clone(),
        });
    }
    if !options.root.is_dir() {
        return Err(ServeError::NoRoot {
            path: options.root.clone(),
        });
    }
    // A project folder itself is not a projects root: its migrations
    // would never list. Say which folder to serve instead.
    if options.root.join(popkern_project::MANIFEST).is_file() {
        return Err(ServeError::ProjectFolder {
            path: options.root.clone(),
        });
    }
    // Absolute, so the page and MCP clients see where projects live.
    let mut options = options;
    if let Ok(root) = options.root.canonicalize() {
        options.root = root;
    }
    let listener = TcpListener::bind(options.bind).map_err(|e| ServeError::Bind {
        addr: options.bind,
        detail: e.to_string(),
    })?;
    Ok((
        listener,
        ServeHandle {
            context: Arc::new(Context { options }),
        },
    ))
}

/// A bound server, ready to accept connections.
pub struct ServeHandle {
    context: Arc<Context>,
}

impl ServeHandle {
    /// Serve `listener` forever (one thread per connection).
    pub fn run(self, listener: TcpListener) {
        for stream in listener.incoming().flatten() {
            let context = Arc::clone(&self.context);
            std::thread::spawn(move || http::connection(stream, &context));
        }
    }
}
