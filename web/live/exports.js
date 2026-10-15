// Downloads of the live model: STEP (AP242, the product carries the
// document's name) and the `.pkn` document, which also carries the
// visibility the user set in the assembly tree. The bytes come from the
// worker and are saved through a blob URL; nothing is uploaded. The
// code dialog (show-as-code) reads the same live model.

/** Export the worker's model as `kind` ("step" or "pkn") named `name`.
 *  For a `.pkn` of an assembly, the tree's `hidden` paths are pushed to
 *  the worker first (the live model was built with the script's own
 *  visibility, `envelope` the one it served). Resolves to the reply. */
export async function exportDocument(host, kind, name, envelope, hidden) {
  if (kind === "pkn" && envelope?.assembly) {
    for (const instance of envelope.assembly.instances) {
      const isHidden = hidden.has(instance.path);
      const scripted = envelope.assembly.instances.find((i) => i.path === instance.path);
      if (isHidden || !scripted.visible) await host.setVisible(instance.path, !isHidden);
    }
  }
  return kind === "step" ? host.exportStep(name) : host.exportPkn();
}

/** Hand `bytes` to the browser as a download named `fileName`. */
export function saveBytes(bytes, fileName, mime = "application/octet-stream") {
  const blob = new Blob([bytes], { type: mime });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = fileName;
  document.body.append(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 2_000);
}

/** The download buttons: build the live model if a pre-built one is on
 *  screen, export it as `kind`, and save it named after the document
 *  ("Planetary gearbox" downloads as planetary_gearbox.step; the STEP
 *  product carries the name too). A failed export shows its outcome. */
export async function downloadModel(page, kind) {
  if (!(await page.ensureLive())) return;
  const name = page.tab.stem;
  const reply = await exportDocument(page.host, kind, name, page.lastGood?.envelope, page.assembly.hidden);
  if (!reply.ok) {
    page.showOutcome(reply.outcome);
    return;
  }
  saveBytes(reply.bytes, `${name}.${kind}`);
}

/** Show-as-code: the live document as canonical Rhai in the code
 *  dialog ("use as script" puts it in the editor). */
export async function showModelCode(page) {
  if (!(await page.ensureLive())) return;
  const reply = await page.host.toRhai();
  if (!reply.ok) {
    page.showOutcome(reply.outcome);
    return;
  }
  document.getElementById("code-text").textContent = reply.code;
  document.getElementById("code-dialog").hidden = false;
}
