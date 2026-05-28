# Pipeline diagrams

Two diagrams that visualize the Phase A pipeline at different levels of detail.

| File | What it shows | When it's useful |
|---|---|---|
| [`01_high_level_flowchart.md`](01_high_level_flowchart.md) | One-page overview: 9 sections, cell counts at each stage, checkpoint files | First read, presentations |
| [`02_detailed_flowchart.md`](02_detailed_flowchart.md) | Function-level call graph, AnnData state evolution, cross-section dependencies | When I want to remember what each step adds to `adata` |

Both are written in Mermaid, which renders natively on GitHub when you click into a `.md` file. It also works in VS Code (with the Mermaid extension) and in JupyterLab.

## Exporting as PNG

If I need static images for a slide deck or a paper figure:

```bash
npm install -g @mermaid-js/mermaid-cli
mmdc -i 01_high_level_flowchart.md -o 01_high_level_flowchart.png --width 2400
mmdc -i 02_detailed_flowchart.md   -o 02_detailed_flowchart.png   --width 2400
```

## Editing

Just edit the markdown. Each ` ```mermaid ` block is the diagram source. Mermaid syntax reference: https://mermaid.js.org/intro/
