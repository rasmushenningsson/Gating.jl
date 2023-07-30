# Gating.jl

Experimental GUI package for gating in Julia.
The main idea is to provide an abstract interface that can be implemented by data providers in order to support:

* Flow cytometry data
* Single cell data

Prototyping is done with Single Cell data using [SingleProjections.jl](https://github.com/rasmushenningsson/SingleCellProjections.jl).


## Usage
Note: Since the package is experimental, everything below can change as the package is updated.

You can downloaded the files used below from [here](https://github.com/rasmushenningsson/SingleCellExampleData).

### Loading data and launching the User Interface

First we load some packages and the Single Cell count data.
```julia
using Gating
using SingleCellProjections
using CSV
using DataFrames

base_path = "/path/to/downloaded/files" # update this path to find the files

# Load count data
sample_path = joinpath(base_path, "GSE164378_RNA_ADT_3P_P1.h5")
counts = load10x(sample_path)
```

Let's annotate the data a bit (this part can be skipped):
```julia
# Annotate each cell with fraction of reads in MT genes
var_counts_fraction!(counts, "name"=>contains(r"^MT-"), "feature_type"=>isequal("Gene Expression"), "fraction_mt")

# Load some more cell annotations
cell_annotations = CSV.read(joinpath(base_path, "GSE164378_RNA_ADT_3P.csv.gz"), DataFrame)
leftjoin!(counts.obs, cell_annotations; on=:barcode)
```

Then we transform and normalize the data
```julia
transformed = sctransform(counts)
normalized = normalize_matrix(transformed)
```

Now we can launch the interactive user interface and select two genes:
```julia
gater = gate(normalized)
gater.x[] = "CD34"
gater.y[] = "CD4"
```

### Interacting directly with the plot
Like any Makie plot, you can pan and zoom to explore the current visualization.
In particular, ctrl+left mouse button will reset the view to show include all the current data.

If you hold down the left shift and click with the left mouse button, you can create a polygon in the plot that is used for selecting cells.
Press escape to remove the current polygon.

After you have created a polygon, press enter to perform gating - i.e. selecting all the cells within the polygon and remove all other cells.
By pressing space, you instead select everyting outside the polygon.


### Using the REPL to manipulate the plot and data
You can also interact with the plot through the REPL.

#### Undo gating
Call `pop!(gater)` to undo the last gating operation.

#### Change axes
To update the x (or y) axes to show another gene, just call e.g.:
```julia
gater.x[] = "TCF3"
```

It is also possible to show principal components (currently supports PCs 1-50):
```julia
gater.x[] = "pc1"
gater.y[] = "pc2" # or "pc3" etc.
```
or a 2d UMAP:
```julia
gater.x[] = "umap1"
gater.y[] = "umap2"
```

Finally, we can also use numerical cell annotations, e.g.:
```julia
gater.x[] = "obs.fraction_mt"
```

#### Coloring
Everything used for the x/axes can also be used for coloring:
```julia
gater.color[] = "CD4" # or "pc1" or "obs.fraction_mt"
```

It is also possible to use categorical annotations for coloring:
```julia
gater.color[] = "celltype.l1"
```

#### Annotations
The current polygon can also be used to create/update cell annotations.
This will update the `DataMatrix` object that was used to launch the user interface.

To create/update a column "MyColumn", with values `true` for cells within the polygon and `false` for those outside, use:
```julia
annotate!(gater, "MyColumn", true, false)
```
Note that cells not included in the current plot are not updated (or set to missing if the column was created.)


Skipping the last argument will only set the values inside the polygon:
```julia
annotate!(gater, "MyColumn", true)
```
