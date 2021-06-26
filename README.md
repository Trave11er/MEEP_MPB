# Scripts for running MPB and MEEP [ARCHIVED]

**DISCLAIMER** This is bad code kept only for reproducibility reasons (see below)

- Scripts to generate .ctl input files for MPB and MEEP electromagnetic codes to simulate photonic topological insulator cavities as in:

   - G. Siroki, P. A. Huidobro and V. Giannini, Phys. Rev. B(Rapid) 96, 041408 (2017)

- based on the photonic topological insulator design described here:

   - L. Wu and X. Hu, Phys. Rev. Lett. 114, 223901 (2015)


## Getting Started

To compile type
```
g++ -o scriptName scriptName.cc
```

### Dependencies

Jmol (optional) can be used to visualize the .cell file representing the structure in the generated .ctl file.

## Authors

* **Gleb Siroki**

## License

This project is licensed under the MIT License - see the LICENSE.md file for details
