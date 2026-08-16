# Native GUI Preprocessor

The native GUI is a two-dimensional RustyFEM model viewer, editor, and early
postprocessor. It is a separate workspace package and uses the public
`rusty-fem` API in the same way as another library consumer.

Run it with the default example:

```bash
cargo run -p rusty-fem-gui
```

Open a specific model at startup:

```bash
cargo run -p rusty-fem-gui -- examples/t3_cantilever.json
```

## Model Workflow

The left panel contains the complete first editing workflow:

- `Open...` opens a JSON model using the native file dialog.
- `Clear canvas` creates an empty model with one material and compatible
  plane-stress, truss, and beam sections.
- `Save` writes to the current path and `Save as...` chooses a new JSON path.
- the `unsaved` status remains visible until the current model is saved.
- `Undo`, `Redo`, and `Delete` edit complete validated model snapshots.
- `Fit view`, pan, zoom, grid visibility, and grid snapping control the canvas.

The JSON written by the GUI is a normal `Model2DInput` file and can be passed to
the command-line solver.

## Drawing FEM Models

Select `Draw FEM` and choose a tool:

- `Insert node`: click the canvas to add an automatically numbered node at the
  displayed coordinates.
- `Insert element`: choose Truss, Beam, T3, T6, Q4, or Q8, then click its corner
  points around the element boundary. The editor creates the element and its
  numbered nodes directly; points coincident with existing nodes reuse them.
  T6 needs three corner clicks and Q8 needs four. Their midside nodes are
  generated automatically and shared by adjacent elements.
- `Move node`: drag an existing node. The resulting geometry is accepted only
  if the complete model remains valid.
- `Constraint`: choose constrained degrees of freedom and click a node. The
  selection replaces that node's previous constraints; unavailable degrees of
  freedom are disabled.
- `Nodal load`: choose an available degree of freedom and value, then click a
  node.
- `Edge traction`: enter global `tx` and `ty`, then click a plane-stress edge.

Creating one element, including all of its new nodes, is one Undo/Redo action.
Middle-button dragging pans while drawing. The mouse wheel zooms. `Escape`
cancels a partially drawn element.

Keyboard commands:

- `Cmd/Ctrl+Z`: Undo.
- `Cmd/Ctrl+Shift+Z` or `Cmd/Ctrl+Y`: Redo.
- `Cmd/Ctrl+S`: Save.
- `Cmd/Ctrl+Shift+S`: Save As.
- `Delete` or `Backspace`: delete the selected node or element.

A node referenced by an element cannot be deleted. Delete the connected element
first. Loads assigned to a deleted element and boundary data assigned to a
deleted free node are removed with their owner.

## Analysis And Results

The `Analysis` section selects the dense or sparse solver and runs the current
model. The calculation runs on a background thread, so the model viewer remains
responsive. During a sparse solve the panel shows the current stage, elapsed
time, CG iteration, absolute and relative residual, and target tolerance. The
progress bar measures logarithmic convergence from a relative residual of one
to the requested tolerance. Result recovery reports the number of completed
elements and advances through the final part of the progress bar. A successful
solve enables these canvas views:

- `Model`: undeformed model and preprocessing symbols.
- `Displacement`: undeformed model plus a scaled deformed outline.
- `sigma x`, `sigma y`, and `tau xy`: signed in-plane stress components.
- `von Mises stress`: equivalent plane-stress intensity.
- `epsilon x`, `epsilon y`, and `gamma xy`: signed in-plane strain components.
- `equivalent strain`: elastic equivalent strain calculated as von Mises stress
  divided by the element material's Young modulus before nodal averaging.

Plane-stress values are recovered at element nodes and arithmetically averaged
between adjacent elements, similarly to common engineering postprocessors. Q4
corner values are extrapolated from Gauss points. T6 values are sampled at all
three corner and three midside nodes, then quadratically interpolated over small
canvas triangles to show variation inside the element.

Scalar fields use the same blue-cyan-green-yellow-orange-red palette as the SVG
renderer. The default legend range follows the current result automatically; it
can be switched to an explicit minimum and maximum for comparing multiple
models on one scale. Selecting a node or plane-stress element reports all
available stress and strain components. The element readout is the mean of its
recovered nodal values. Sparse runs also show iteration count, relative residual,
and termination status. Beam and truss elements currently show deformation but
do not receive a plane-stress contour.

Any model edit invalidates the current results, so stale contours cannot remain
visible after geometry, boundary conditions, or loads change.

## Module Layout

The repository is a Cargo workspace with two packages:

- the root `rusty-fem` package contains the reusable FEM library and CLI;
- `crates/rusty-fem-gui` contains the native application and depends on the
  root package through a local path dependency.

The GUI implementation lives under `crates/rusty-fem-gui/src/gui`:

- `app.rs`: application state, panels, and validated edit commands.
- `canvas.rs`: projection, drawing, picking, direct manipulation, and result layers.
- `document.rs`: model snapshots, Undo/Redo, and JSON persistence.
- `results.rs`: solver execution and compact postprocessing data.
- `loaded_model.rs`: model loading and coordinate bounds.
- `model_browser.rs`: discovery of example JSON files.
- `selection.rs`: selected entities and view toggles.
- `topology.rs`: visible outlines and higher-order edge segmentation.
- `workflow.rs`: work modes, drawing tools, and element draft types.
- `theme.rs`: shared native GUI styling.

## Data Flow

```text
JSON -> Model2DInput -> validated Model2D -> canvas
                         |       ^
                         v       |
                    edit snapshot + validation -> Undo/Redo
                         |
                         +-> JSON Save
                         +-> dense/sparse solve -> displacement/stress views
```

`Model2D` remains the canonical solver model. GUI edits operate on a serializable
snapshot and rebuild the model through `Model2DInput::into_model`, so duplicate
IDs, invalid connectivity, incompatible sections, and degenerate geometry are
rejected before an edit becomes part of the document history. Drawing assigns
the next available element and node IDs automatically and reuses coincident
corner or midside nodes within the editor tolerance.
