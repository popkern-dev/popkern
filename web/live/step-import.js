// Open STEP (PK-631): a STEP file opens as a view-only assembly — the
// product tree, names and placements, with no script behind it. The
// file is read in the page and handed to the worker; it never leaves
// the tab.

/** Open `bytes` (the file `fileName`) on `page`. Resolves to the
 *  worker's reply, or null while the kernel is still loading. */
export async function openStepAssembly(page, fileName, bytes) {
  if (!page.host.info) {
    page.status("the kernel is still loading", true);
    return null;
  }
  await page.runDone;
  const stem = fileName.replace(/\.(step|stp)$/i, "");
  page.status(`opening ${fileName}…`);
  const token = ++page.viewToken;
  page.certifyGeneration += 1;
  const reply = await page.host.openStep(bytes, stem);
  if (token !== page.viewToken) return reply;
  if (reply.ok && reply.envelope?.outcome === "served") {
    // The editor shows a placeholder the run path recognises (and
    // declines to build) until the user picks or writes a script.
    page.tab.importedSource = `// ${fileName}: an imported STEP assembly — view only (no script).\n// Pick an example or write a script to build.\n`;
    page.tab.open(stem, page.tab.importedSource);
    page.viewer.fitted = false;
    page.live = true;
  }
  await page.receive(reply);
  if (reply.ok && reply.envelope?.printed?.length) page.status(`${reply.envelope.printed[0]} — ${page.servedText ?? ""}`);
  return reply;
}
