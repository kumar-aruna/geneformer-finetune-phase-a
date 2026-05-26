# Pipeline Diagrams

Two complementary visualizations of the Phase A scRNA-seq analysis pipeline.

| File | What it shows | When to use |
|---|---|---|
| [`01_high_level_flowchart.md`](01_high_level_flowchart.md) | One-page overview — 9 sections, cell counts, checkpoint files | First read; presentations; high-level discussion |
| [`02_detailed_flowchart.md`](02_detailed_flowchart.md) | Function-level call graph + AnnData state evolution + cross-section dependencies | Interview prep; understanding what each step does to `adata` |

Both use **Mermaid** — they render natively in:
- GitHub README and markdown files
- VS Code (with the Mermaid extension)
- JupyterLab (via `mermaid` cell magic or `IPython.display`)
- Any modern markdown viewer

## To export as PNG

If you need static images (e.g., for a slide deck or paper figure):

```bash
# Install the Mermaid CLI once:
npm install -g @mermaid-js/mermaid-cli

# Then export both diagrams:
mmdc -i 01_high_level_flowchart.md -o 01_high_level_flowchart.png --width 2400
mmdc -i 02_detailed_flowchart.md   -o 02_detailed_flowchart.png   --width 2400
```

## To edit

Edit the markdown directly. Each `\`\`\`mermaid` block contains a diagram source. Mermaid syntax is text-based — easy to version-control and review changes through pull requests.

For the Mermaid syntax reference: https://mermaid.js.org/intro/
