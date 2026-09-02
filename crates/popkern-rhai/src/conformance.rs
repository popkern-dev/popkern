//! The front-end conformance checklist (design §13): a machine-readable
//! inventory generated from `docs/architecture/scripting-surface.md` —
//! every op with its stability mark and parameter names, the selection
//! predicates, and the §10 error classes with their payload fields.
//! Skins are audited against it (`conformance_audit`), so a spec row
//! without an implementation fails CI instead of drifting.

use serde::{Deserialize, Serialize};

/// One operation row of the spec.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct OpRow {
    /// Operation name (`fillet`).
    pub name: String,
    /// Spec section (`5`, `4.1`, `7`).
    pub section: String,
    /// `stable`, `preview`, or `reserved`.
    pub mark: String,
    /// The canonical signature text, when the row gives one.
    pub signature: String,
    /// Parameter names in signature order.
    pub parameters: Vec<Parameter>,
}

/// A signature parameter.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Parameter {
    /// Name as written.
    pub name: String,
    /// `true` when spelled with a trailing `?`.
    pub optional: bool,
}

/// A selection predicate of the §6.1 grammar.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Predicate {
    /// Predicate name (`between`).
    pub name: String,
    /// `edge` or `face`.
    pub kind: String,
    /// `stable` or `preview`.
    pub mark: String,
}

/// One §10 error class.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ErrorClass {
    /// Class name (`BlendDeclined`).
    pub name: String,
    /// Spec section (`10.2` / `10.3`).
    pub section: String,
    /// Payload field names the row lists.
    pub fields: Vec<String>,
}

/// The generated checklist.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Checklist {
    /// The spec the checklist was generated from.
    pub source: String,
    /// Every operation row, in document order.
    pub ops: Vec<OpRow>,
    /// The selection grammar.
    pub predicates: Vec<Predicate>,
    /// The error taxonomy.
    pub error_classes: Vec<ErrorClass>,
}

/// Predicates the §6.1 prose marks `preview` (the four the recorder
/// added over the `popkern-part` vocabulary).
const PREVIEW_PREDICATES: [&str; 4] = ["conic", "nurbs", "labeled", "named"];

fn backticked(cell: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut rest = cell;
    while let Some(start) = rest.find('`') {
        let after = &rest[start + 1..];
        let Some(end) = after.find('`') else { break };
        out.push(after[..end].to_string());
        rest = &after[end + 1..];
    }
    out
}

fn is_identifier(text: &str) -> bool {
    !text.is_empty()
        && text
            .chars()
            .all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || c == '_')
}

/// Split a table row into cells, honouring `\|` escapes.
fn cells(line: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut current = String::new();
    let mut chars = line.trim().trim_start_matches('|').chars().peekable();
    while let Some(c) = chars.next() {
        match c {
            '\\' if chars.peek() == Some(&'|') => {
                current.push('|');
                chars.next();
            }
            '|' => {
                out.push(current.trim().to_string());
                current.clear();
            }
            other => current.push(other),
        }
    }
    if !current.trim().is_empty() {
        out.push(current.trim().to_string());
    }
    out
}

/// Parameter names of `name(a, b?, c)` from a signature text.
fn parameters_of(signature: &str, name: &str) -> Option<Vec<Parameter>> {
    let needle = format!("{name}(");
    let start = signature.find(&needle)? + needle.len();
    let mut depth = 0usize;
    let mut end = None;
    for (i, c) in signature[start..].char_indices() {
        match c {
            '(' | '[' | '{' => depth += 1,
            ')' | ']' | '}' if depth == 0 => {
                end = Some(start + i);
                break;
            }
            ')' | ']' | '}' => depth -= 1,
            _ => {}
        }
    }
    let inner = &signature[start..end?];
    let mut params = Vec::new();
    let mut depth = 0usize;
    let mut current = String::new();
    for c in inner.chars() {
        match c {
            '(' | '[' | '{' => {
                depth += 1;
                current.push(c);
            }
            ')' | ']' | '}' => {
                depth -= 1;
                current.push(c);
            }
            ',' if depth == 0 => {
                params.push(std::mem::take(&mut current));
            }
            _ => current.push(c),
        }
    }
    if !current.trim().is_empty() {
        params.push(current);
    }
    Some(
        params
            .iter()
            .filter_map(|p| {
                let p = p.trim();
                let optional = p.ends_with('?');
                let head = p.trim_end_matches('?');
                let name = head.split([':', '=', ' ']).next()?.trim();
                is_identifier(name).then(|| Parameter {
                    name: name.to_string(),
                    optional,
                })
            })
            .collect(),
    )
}

struct Table {
    section: String,
    header: Vec<String>,
    rows: Vec<Vec<String>>,
}

fn tables(markdown: &str) -> (Vec<Table>, Vec<(String, Vec<String>)>) {
    let mut section = String::new();
    let mut out = Vec::new();
    let mut current: Option<Table> = None;
    let mut blocks: Vec<(String, Vec<String>)> = Vec::new();
    let mut in_block: Option<Vec<String>> = None;
    for line in markdown.lines() {
        if let Some(rest) = line
            .strip_prefix("### ")
            .or_else(|| line.strip_prefix("## "))
        {
            section = rest
                .split_whitespace()
                .next()
                .unwrap_or("")
                .trim_end_matches('.')
                .to_string();
        }
        if line.trim_start().starts_with("```") {
            match in_block.take() {
                Some(lines) => blocks.push((section.clone(), lines)),
                None => in_block = Some(Vec::new()),
            }
            continue;
        }
        if let Some(lines) = in_block.as_mut() {
            lines.push(line.to_string());
            continue;
        }
        if line.trim_start().starts_with('|') {
            let row = cells(line);
            if row
                .iter()
                .all(|c| c.chars().all(|ch| ch == '-' || ch == ':'))
            {
                continue;
            }
            match current.as_mut() {
                Some(table) => table.rows.push(row),
                None => {
                    current = Some(Table {
                        section: section.clone(),
                        header: row,
                        rows: Vec::new(),
                    })
                }
            }
        } else if let Some(table) = current.take() {
            out.push(table);
        }
    }
    if let Some(table) = current.take() {
        out.push(table);
    }
    (out, blocks)
}

fn op_rows(table: &Table, default_mark: &str) -> Vec<OpRow> {
    let column = |needle: &str| {
        table
            .header
            .iter()
            .position(|h| h.to_ascii_lowercase().contains(needle))
    };
    let mark_col = column("mark");
    let signature_col = column("signature");
    let mut rows = Vec::new();
    for row in &table.rows {
        let Some(first) = row.first() else { continue };
        if first.starts_with('↳') {
            continue;
        }
        let names: Vec<String> = backticked(first)
            .into_iter()
            .filter(|n| is_identifier(n))
            .collect();
        if names.is_empty() {
            continue;
        }
        let marks: Vec<String> = mark_col
            .and_then(|i| row.get(i))
            .map(|m| m.split('/').map(|s| s.trim().to_string()).collect())
            .unwrap_or_else(|| vec![default_mark.to_string()]);
        let signature_cell = signature_col
            .and_then(|i| row.get(i))
            .cloned()
            .unwrap_or_default();
        for (i, name) in names.iter().enumerate() {
            let mark = marks
                .get(i)
                .or_else(|| marks.first())
                .cloned()
                .unwrap_or_else(|| default_mark.to_string());
            let signature = backticked(&signature_cell)
                .into_iter()
                .find(|s| s.contains(&format!("{name}(")))
                .unwrap_or_default();
            let parameters = parameters_of(&signature, name).unwrap_or_default();
            rows.push(OpRow {
                name: name.clone(),
                section: table.section.clone(),
                mark,
                signature,
                parameters,
            });
        }
    }
    rows
}

fn error_rows(table: &Table) -> Vec<ErrorClass> {
    table
        .rows
        .iter()
        .filter_map(|row| {
            let first = row.first()?;
            if first.starts_with('↳') {
                return None;
            }
            let name = backticked(first).into_iter().next()?;
            let payload = row.get(1).cloned().unwrap_or_default();
            let fields: Vec<String> = backticked(&payload)
                .into_iter()
                .map(|f| f.trim_end_matches("[]").to_string())
                .filter(|f| is_identifier(f))
                .collect();
            Some(ErrorClass {
                name,
                section: table.section.clone(),
                fields,
            })
        })
        .collect()
}

/// Builder-style ops of a code block: `.name(` occurrences, plus a
/// `model.name(` entry op.
fn block_ops(section: &str, lines: &[String], mark: &str) -> Vec<OpRow> {
    let mut seen = std::collections::BTreeSet::new();
    let mut rows = Vec::new();
    for line in lines {
        let mut rest = line.as_str();
        while let Some(dot) = rest.find('.') {
            let after = &rest[dot + 1..];
            let name: String = after
                .chars()
                .take_while(|c| c.is_ascii_lowercase() || *c == '_')
                .collect();
            if !name.is_empty() && after[name.len()..].starts_with('(') && seen.insert(name.clone())
            {
                let signature =
                    after[..after.find(')').map(|i| i + 1).unwrap_or(after.len())].to_string();
                rows.push(OpRow {
                    parameters: parameters_of(&signature, &name).unwrap_or_default(),
                    name,
                    section: section.to_string(),
                    mark: mark.to_string(),
                    signature,
                });
            }
            rest = after;
        }
    }
    rows
}

fn grammar_predicates(lines: &[String]) -> Vec<Predicate> {
    let mut out = Vec::new();
    let mut kind = "edge";
    for line in lines {
        if line.contains("face predicate") {
            kind = "face";
        }
        if line.contains(":=") && !line.contains("predicate") {
            continue;
        }
        let mut rest = line.as_str();
        while let Some(open) = rest.find('(') {
            let head: String = rest[..open]
                .chars()
                .rev()
                .take_while(|c| c.is_ascii_lowercase() || *c == '_')
                .collect::<Vec<_>>()
                .into_iter()
                .rev()
                .collect();
            if !head.is_empty() && !matches!(head.as_str(), "edges" | "faces" | "or") {
                let mark = if PREVIEW_PREDICATES.contains(&head.as_str()) {
                    "preview"
                } else {
                    "stable"
                };
                if !out
                    .iter()
                    .any(|p: &Predicate| p.name == head && p.kind == kind)
                {
                    out.push(Predicate {
                        name: head,
                        kind: kind.to_string(),
                        mark: mark.to_string(),
                    });
                }
            }
            rest = &rest[open + 1..];
        }
    }
    out
}

fn prose_ops(markdown: &str, section: &str, mark: &str) -> Vec<OpRow> {
    let heading = format!("### {section} ");
    let Some(start) = markdown.find(&heading) else {
        return Vec::new();
    };
    let body = &markdown[start..];
    let end = body[heading.len()..]
        .find("\n## ")
        .or_else(|| body[heading.len()..].find("\n### "))
        .map(|i| i + heading.len())
        .unwrap_or(body.len());
    let mut seen = std::collections::BTreeSet::new();
    backticked(&body[..end])
        .into_iter()
        .filter_map(|token| {
            let name = token.split('(').next()?.trim().to_string();
            (is_identifier(&name) && seen.insert(name.clone())).then(|| OpRow {
                parameters: parameters_of(&token, &name).unwrap_or_default(),
                name,
                section: section.to_string(),
                mark: mark.to_string(),
                signature: token.clone(),
            })
        })
        .collect()
}

/// Generate the checklist from the spec's markdown.
pub fn checklist_from_spec(markdown: &str, source: &str) -> Checklist {
    let (tables, blocks) = tables(markdown);
    let mut ops = Vec::new();
    let mut error_classes = Vec::new();
    let mut predicates = Vec::new();
    for (section, lines) in &blocks {
        match section.as_str() {
            "3.1" => ops.extend(block_ops("3.1", lines, "stable")),
            "3.2" => ops.extend(block_ops("3.2", lines, "preview")),
            "6.1" => predicates = grammar_predicates(lines),
            _ => {}
        }
    }
    ops.extend(prose_ops(markdown, "4.3", "preview"));
    for table in &tables {
        let first = table.header.first().map(String::as_str).unwrap_or("");
        match first {
            "Op" | "Verb" => {
                let default = if table.section == "4.1" {
                    "stable"
                } else {
                    "preview"
                };
                ops.extend(op_rows(table, default));
            }
            "Class" => error_classes.extend(error_rows(table)),
            _ => {}
        }
    }
    ops.sort_by(|a, b| a.section.cmp(&b.section).then(a.name.cmp(&b.name)));
    ops.dedup_by(|a, b| a.name == b.name && a.section == b.section);
    Checklist {
        source: source.to_string(),
        ops,
        predicates,
        error_classes,
    }
}
