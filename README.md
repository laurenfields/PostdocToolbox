# PostdocToolbox

Personal toolbox of small utilities for proteomics and data workflows.
A single home for the one-off Python scripts I keep writing for DIA-NN /
Skyline / mass-spec pipelines so they stop scattering across project
directories.

## Layout

Tools are grouped by domain. Each tool gets its own subdirectory with its
own short README and script(s).

```
PostdocToolbox/
├── README.md
├── LICENSE
├── .gitignore
└── proteomics/
    └── diann_lib_to_openswath_tsv/
        ├── README.md
        └── convert_to_openswath.py
```

Domain folders that may show up over time: `proteomics/`,
`bioinformatics/`, `data_wrangling/`, etc.

## Adding a new tool

Pick (or create) the appropriate domain folder, add a subdirectory named
after the tool, and drop a `README.md` plus the script(s) inside. The
README should cover what / why / dependencies / usage and any version
caveats. Keep each tool self-contained and standalone — these are not part
of a Python package.

## License

MIT. See [LICENSE](LICENSE).
