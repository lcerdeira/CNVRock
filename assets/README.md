# assets

Small static files that can be committed to the repository. Typical contents:

- Figures and diagrams embedded in documentation or `README.md` files
- Logo / banner images
- Small reference files (e.g. colour palettes, font files)
- Example input/output samples used in notebooks or demos

## Conventions

- Keep file sizes small (< 1 MB per file as a rule of thumb). Large binary
  assets belong in external storage, not the repo.
- Use descriptive filenames: `architecture_diagram.png`, not `fig1.png`.
- Where applicable, store the source file alongside the export
  (e.g. `diagram.drawio` next to `diagram.png`).
