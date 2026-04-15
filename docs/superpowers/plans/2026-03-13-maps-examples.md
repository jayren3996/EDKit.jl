# Maps Examples Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a dedicated `examples/Maps/` notebook sequence that explains `DoubleBasis`, `symmetrizer`, and rectangular operator construction.

**Architecture:** Create a new topic folder with three notebook-first examples and a local README. Each notebook should have explicit explanatory markdown sections, small physically meaningful examples, and runnable code that loads the local package. Then update the top-level example indexes so the new topic is discoverable.

**Tech Stack:** Julia, Jupyter notebooks (`.ipynb`), existing EDKit example-loading pattern, README markdown.

---

## Chunk 1: Notebook Set

### Task 1: Add the folder and index

**Files:**
- Create: `/Users/ren/Library/CloudStorage/OneDrive-UniversityofLeeds/GitHub/EDKit.jl/examples/Maps/README.md`

- [ ] **Step 1: Check that the folder does not already exist**

Run: `find examples -maxdepth 2 -type d | sort | rg 'examples/Maps$'`
Expected: no output

- [ ] **Step 2: Create the folder index**

Add a README that explains the topic split and lists:
- `DoubleBasisBasics.ipynb`
- `Symmetrizers.ipynb`
- `InterbasisOperators.ipynb`

- [ ] **Step 3: Verify the README exists**

Run: `sed -n '1,200p' examples/Maps/README.md`
Expected: markdown index with the three notebooks

### Task 2: Add the notebooks

**Files:**
- Create: `/Users/ren/Library/CloudStorage/OneDrive-UniversityofLeeds/GitHub/EDKit.jl/examples/Maps/DoubleBasisBasics.ipynb`
- Create: `/Users/ren/Library/CloudStorage/OneDrive-UniversityofLeeds/GitHub/EDKit.jl/examples/Maps/Symmetrizers.ipynb`
- Create: `/Users/ren/Library/CloudStorage/OneDrive-UniversityofLeeds/GitHub/EDKit.jl/examples/Maps/InterbasisOperators.ipynb`

- [ ] **Step 1: Write a failing existence check**

Run:
```bash
for f in \
  examples/Maps/DoubleBasisBasics.ipynb \
  examples/Maps/Symmetrizers.ipynb \
  examples/Maps/InterbasisOperators.ipynb; do
  test -f "$f"
done
```
Expected: non-zero exit because the notebooks do not exist yet

- [ ] **Step 2: Create `DoubleBasisBasics.ipynb`**

Include sections:
- what `DoubleBasis(Btarget, Bsource)` stores
- direction of the map
- fast `T(v)` action versus `symmetrizer(T) * v`
- an `AbelianBasis <- TensorBasis` example first

- [ ] **Step 3: Create `Symmetrizers.ipynb`**

Include sections:
- full-space state to sector coordinates
- lifted projected state `P * (P' * ψ)`
- parity-sector example and one Abelian/momentum example

- [ ] **Step 4: Create `InterbasisOperators.ipynb`**

Include sections:
- rectangular `operator(..., B::DoubleBasis)`
- `trans_inv_operator(...)` with `DoubleBasis`
- matrix shapes and physical interpretation

- [ ] **Step 5: Run each notebook in script mode**

Run the local notebook runner against each file and confirm success.

## Chunk 2: Docs Wiring

### Task 3: Update example indexes

**Files:**
- Modify: `/Users/ren/Library/CloudStorage/OneDrive-UniversityofLeeds/GitHub/EDKit.jl/examples/README.md`
- Modify: `/Users/ren/Library/CloudStorage/OneDrive-UniversityofLeeds/GitHub/EDKit.jl/README.md`

- [ ] **Step 1: Add `Maps/` to the examples topic list**

- [ ] **Step 2: Add the three notebooks to the recommended starting points**

- [ ] **Step 3: Verify links and wording**

Run:
```bash
sed -n '1,220p' examples/README.md
sed -n '1,260p' README.md
```
Expected: `Maps/` appears as a first-class topic and the notebook names are listed correctly

## Chunk 3: Validation

### Task 4: Run validation

**Files:**
- Validate only

- [ ] **Step 1: Re-run the three new notebooks**

- [ ] **Step 2: Re-run `Pkg.test()` only if notebook changes require package-code validation**

- [ ] **Step 3: Summarize any remaining coverage gap**

The summary should note whether there is still a missing map-related notebook, especially around MPS projection versus basis-to-basis maps.
