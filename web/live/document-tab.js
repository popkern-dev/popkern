// The open document: its title (which names the editor tab and every
// download), the text it was opened with (so the tab can mark an edit),
// and whether it is an imported STEP assembly with no script behind it.

const $ = (id) => document.getElementById(id);

/** The default script's title (the gallery titles every example). */
export const DEFAULT_TITLE = "Boss plate";
/** The title of a script that is neither the default nor an example. */
export const UNTITLED = "Document";

/** A title as a file stem: lowercase, spaces as underscores, and no
 *  characters a file name cannot carry ("Planetary gearbox" ->
 *  planetary_gearbox). */
export function fileStem(title) {
  const stem = title
    .trim()
    .toLowerCase()
    .replace(/\s+/g, "_")
    .replace(/[\/\\:*?"<>|\x00-\x1f]/g, "");
  return stem || "document";
}

export class DocumentTab {
  /** `editor` is the page's Editor: the tab compares its text with the
   *  text the document was opened with. */
  constructor(editor) {
    this.editor = editor;
    this.title = null;
    this.source = null;
    // The placeholder text of an imported STEP document (view only), or
    // null when the document is a script.
    this.importedSource = null;
  }

  /** Put `text` in the editor as the document titled `title`. */
  open(title, text) {
    this.title = title;
    this.source = text;
    this.editor.value = text;
    this.show();
  }

  /** The tab: the document's name, marked once it is edited. */
  show() {
    const edited = this.source != null && this.editor.value !== this.source;
    $("doc-name").textContent = `${this.stem}.rhai${edited ? " ●" : ""}`;
    $("doc-name").title = edited ? "edited: builds live" : "";
  }

  /** The document's file stem, for the tab and downloads. */
  get stem() {
    return fileStem(this.title ?? UNTITLED);
  }

  /** Whether `source` is the imported STEP document's placeholder. */
  isImported(source) {
    return this.importedSource != null && source === this.importedSource;
  }
}
