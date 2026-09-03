//! PK-413: the playground serves the scripting guide as a page.
//! `web/playground/guide.html` is generated from
//! `docs/reference/scripting-guide.md` by the small Markdown renderer
//! below (headings, paragraphs, lists, tables, fenced code, inline
//! code/bold/links — the subset the guide uses; no dependency).
//! `POPKERN_UPDATE_GALLERY=1` regenerates it; otherwise the test fails
//! on drift so the page and the guide cannot disagree.

use std::path::{Path, PathBuf};

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("workspace root")
}

fn escape(text: &str) -> String {
    text.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

/// Inline markdown: code spans, bold, links.
fn inline(text: &str) -> String {
    let mut out = String::new();
    let mut rest = text;
    while !rest.is_empty() {
        if let Some(start) = rest.find('`') {
            out.push_str(&inline_plain(&rest[..start]));
            let after = &rest[start + 1..];
            if let Some(end) = after.find('`') {
                out.push_str(&format!("<code>{}</code>", escape(&after[..end])));
                rest = &after[end + 1..];
            } else {
                out.push_str(&inline_plain(after));
                break;
            }
        } else {
            out.push_str(&inline_plain(rest));
            break;
        }
    }
    out
}

fn inline_plain(text: &str) -> String {
    let mut out = escape(text);
    // **bold**
    while let Some(start) = out.find("**") {
        let after = &out[start + 2..];
        let Some(end) = after.find("**") else { break };
        let bold = format!("<strong>{}</strong>", &after[..end]);
        out = format!("{}{}{}", &out[..start], bold, &after[end + 2..]);
    }
    // [text](href)
    while let Some(start) = out.find('[') {
        let Some(mid) = out[start..].find("](") else {
            break;
        };
        let Some(end) = out[start + mid..].find(')') else {
            break;
        };
        let text = &out[start + 1..start + mid];
        let href = &out[start + mid + 2..start + mid + end];
        let link = format!("<a href=\"{href}\">{text}</a>");
        out = format!("{}{}{}", &out[..start], link, &out[start + mid + end + 1..]);
    }
    out
}

/// Render the guide's Markdown subset to an HTML body.
fn render(markdown: &str) -> String {
    let mut html = String::new();
    let mut lines = markdown.lines().peekable();
    let mut paragraph: Vec<String> = Vec::new();
    let flush = |paragraph: &mut Vec<String>, html: &mut String| {
        if !paragraph.is_empty() {
            html.push_str(&format!("<p>{}</p>\n", inline(&paragraph.join(" "))));
            paragraph.clear();
        }
    };
    while let Some(line) = lines.next() {
        if let Some(fence) = line.strip_prefix("```") {
            flush(&mut paragraph, &mut html);
            let lang = fence.trim();
            let mut code = String::new();
            for l in lines.by_ref() {
                if l.starts_with("```") {
                    break;
                }
                code.push_str(l);
                code.push('\n');
            }
            html.push_str(&format!(
                "<pre class=\"lang-{}\"><code>{}</code></pre>\n",
                if lang.is_empty() { "text" } else { lang },
                escape(&code)
            ));
        } else if let Some(h) = line.strip_prefix("### ") {
            flush(&mut paragraph, &mut html);
            html.push_str(&format!("<h3>{}</h3>\n", inline(h)));
        } else if let Some(h) = line.strip_prefix("## ") {
            flush(&mut paragraph, &mut html);
            html.push_str(&format!("<h2>{}</h2>\n", inline(h)));
        } else if let Some(h) = line.strip_prefix("# ") {
            flush(&mut paragraph, &mut html);
            html.push_str(&format!("<h1>{}</h1>\n", inline(h)));
        } else if line.starts_with('|') {
            flush(&mut paragraph, &mut html);
            let mut rows: Vec<&str> = vec![line];
            while let Some(next) = lines.peek() {
                if next.starts_with('|') {
                    rows.push(lines.next().unwrap_or_default());
                } else {
                    break;
                }
            }
            html.push_str("<table>\n");
            for (i, row) in rows.iter().enumerate() {
                if i == 1 && row.contains("---") {
                    continue;
                }
                let tag = if i == 0 { "th" } else { "td" };
                let cells: Vec<&str> = row.trim_matches('|').split('|').collect();
                html.push_str("<tr>");
                for cell in cells {
                    html.push_str(&format!("<{tag}>{}</{tag}>", inline(cell.trim())));
                }
                html.push_str("</tr>\n");
            }
            html.push_str("</table>\n");
        } else if line.starts_with("- ")
            || line.starts_with("1. ")
            || line.starts_with("2. ")
            || line.starts_with("3. ")
            || line.starts_with("4. ")
            || line.starts_with("5. ")
        {
            flush(&mut paragraph, &mut html);
            let ordered = !line.starts_with("- ");
            let tag = if ordered { "ol" } else { "ul" };
            let mut items: Vec<String> = Vec::new();
            let mut current = line
                .split_once(' ')
                .map(|x| x.1)
                .unwrap_or_default()
                .to_string();
            while let Some(next) = lines.peek() {
                if next.starts_with("  ") && !next.trim().is_empty() {
                    current.push(' ');
                    current.push_str(next.trim());
                    lines.next();
                } else if next.starts_with("- ")
                    || next.chars().next().is_some_and(|c| c.is_ascii_digit())
                        && next.contains(". ")
                {
                    items.push(current.clone());
                    current = next
                        .split_once(' ')
                        .map(|x| x.1)
                        .unwrap_or_default()
                        .to_string();
                    lines.next();
                } else {
                    break;
                }
            }
            items.push(current);
            html.push_str(&format!("<{tag}>\n"));
            for item in items {
                html.push_str(&format!("<li>{}</li>\n", inline(&item)));
            }
            html.push_str(&format!("</{tag}>\n"));
        } else if line.trim().is_empty() {
            flush(&mut paragraph, &mut html);
        } else {
            paragraph.push(line.trim().to_string());
        }
    }
    flush(&mut paragraph, &mut html);
    html
}

fn page(body: &str) -> String {
    format!(
        "<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n<meta http-equiv=\"Content-Security-Policy\" content=\"default-src 'self'; script-src 'none'; style-src 'self'; img-src 'self'; object-src 'none'; base-uri 'none'; form-action 'none'\">\n<meta name=\"referrer\" content=\"no-referrer\">\n<title>popkern scripting guide</title>\n<link rel=\"icon\" href=\"./favicon.svg\" type=\"image/svg+xml\">\n<link rel=\"stylesheet\" href=\"./guide.css\">\n</head>\n<body>\n<nav class=\"top\"><a href=\"./index.html\">← playground</a></nav>\n<main class=\"guide\">\n{body}</main>\n</body>\n</html>\n"
    )
}

#[test]
fn the_guide_page_is_the_rendered_guide() {
    let markdown =
        std::fs::read_to_string(root().join("docs/reference/scripting-guide.md")).expect("guide");
    let html = page(&render(&markdown));
    let path = root().join("web/playground/guide.html");
    if std::env::var("POPKERN_UPDATE_GALLERY").is_ok_and(|v| v == "1") {
        std::fs::write(&path, &html).expect("write guide.html");
        return;
    }
    let committed = std::fs::read_to_string(&path).unwrap_or_default();
    assert_eq!(
        committed, html,
        "web/playground/guide.html drifted from docs/reference/scripting-guide.md: run with \
         POPKERN_UPDATE_GALLERY=1"
    );
    // Sanity of the render: every heading and every code block made it.
    assert_eq!(
        markdown.lines().filter(|l| l.starts_with("## ")).count(),
        html.matches("<h2>").count()
    );
    assert_eq!(
        markdown
            .lines()
            .filter(|l| l.starts_with("```rhai"))
            .count(),
        html.matches("<pre class=\"lang-rhai\">").count()
    );
}
