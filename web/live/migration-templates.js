// Migration templates (no DOM): the same text as popkern-project's
// `scaffold` module, so a migration the page adds is the one
// `popkern` tools would write (tests/migration_templates.rs compares
// them under node).

/** The first migration of a new project called `name`. */
export function firstMigration(name) {
  return (
    `// ${name}: the first migration builds the document. Later\n` +
    "// migrations (002_<name>.rhai, …) start with `let m = previous();`\n" +
    "// and edit it.\n" +
    "let m = model();\n" +
    'm.profile(xy()).rect(0.08, 0.05).sketch("plate").extrude(blind(0.01));\n' +
    "m\n"
  );
}

/** The source of a new later migration. */
export const NEXT_MIGRATION = "// Edits the document the earlier migrations built.\nlet m = previous();\n\nm\n";

/** A Rhai string literal for `text`. */
export function rhaiString(text) {
  let out = '"';
  for (const c of text) {
    if (c === "\\") out += "\\\\";
    else if (c === '"') out += '\\"';
    else if (c === "\n") out += "\\n";
    else out += c;
  }
  return `${out}"`;
}

/** "Insert a feature here": a migration inserting one feature right
 *  after the step labelled `label` (by label when it is `unique` in
 *  the tree, else by its 1-based `position`). */
export function insertAfterMigration(label, position, unique) {
  const anchor = unique ? `m.feature(${rhaiString(label)})` : `m.features()[${Math.max(0, position - 1)}]`;
  return (
    `// Inserts one feature after ${rhaiString(label)} (step ${position} of the history).\n` +
    "// Record exactly one feature on `body`, for example\n" +
    "//   body.chamfer(edges().at_height(1.0), 0.05)\n" +
    "// The features after it replay after the new one.\n" +
    "let m = previous();\n" +
    `m.insert_after(${anchor}, |m, body| {\n    body\n});\n` +
    "m\n"
  );
}

/** The file stem for an insert-after migration. */
export function insertAfterName(label) {
  let cleaned = "";
  for (const c of label.toLowerCase()) {
    const d = /\s/.test(c) ? "_" : c;
    if (/^[a-z0-9_-]$/.test(d)) cleaned += d;
  }
  return cleaned ? `insert_after_${cleaned}` : "insert_feature";
}
