using Documenter, Optikos

#makedocs(sitename="My Documentation", remotes = nothing)

makedocs(sitename="Optikos.jl")
deploydocs(
    repo = "github.com/lsga001/Optikos.jl.git",
)

pages = [
  "Home" => "index.md",
  "Getting started" => "example.md",
]
