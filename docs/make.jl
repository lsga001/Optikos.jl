using Documenter, Optikos

DocMeta.setdocmeta!(
    Optikos,
    :DocTestSetup,
    :(using Optikos);
    recursive = true
)

makedocs(
  sitename = "Optikos.jl", 
  modules = [Optikos],
  checkdocs = :exports,
  pages = [
    "Home" => "index.md",
    "Examples" => "examples.md",
    "Types" => "types.md",
    "Beams" => "beams.md",
    "Ensembles" => "ensembles.md",
    "Operations" => "operations.md",
    "Optical Elements" => "opticalelements.md",
    "Visualization" => "visualization.md",
  ],
  repo = Remotes.GitHub("lsga001", "Optikos.jl"),
)
deploydocs(
  repo = "github.com/lsga001/Optikos.jl.git",
  devbranch = "main",  
)
