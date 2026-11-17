# This script installs zsh completions for FlexibleSUSY's configure,
# createmodel, createaddon scripts as well as all spectrum generators.
#
# Usage:
#
#   . utils/install-zsh_completions.zsh
#
# Run it from the FlexibleSUSY directory.  Model and addon names are
# looked up relative to the current directory when completing.

if (( ! $+functions[compdef] )); then
    autoload -Uz compinit && compinit
fi

_flexiblesusy_available_models()
{
    local -a models=( model_files/*/FlexibleSUSY.m.in(N:h:t) )
    compadd -a models
}

_flexiblesusy_created_models()
{
    local -a models=( models/*/module.mk(N:h:t) )
    _values -s , model all $models
}

_flexiblesusy_addons()
{
    local -a addons=( addons/*(N/:t) )
    (( $#addons )) && _values -s , addon $addons
}

# The shell doesn't expand ~ in --option=~/path and FlexibleSUSY doesn't
# either, so complete such a path to one which starts with $HOME.
_flexiblesusy_files()
{
    [[ ${PREFIX} == \~/* ]] && PREFIX=${HOME}${PREFIX#\~}
    _files "$@"
}

_flexiblesusy_configure()
{
    local feature
    local -a features args

    features=(
        'colors[Colored output]'
        'compile[Compile the source code]'
        'compiler-warnings[Enable compiler warnings]'
        'debug[Debug messages and assertions]'
        'fflite[Use fflite to compute the loop functions]'
        'gm2calc[Enable GM2Calc]'
        'higgstools[Enable HiggsTools]'
        'himalaya[Enable Himalaya]'
        'lilith[Enable Lilith]'
        'librarylink[Build LibraryLink library for Mathematica interface]'
        'looptools[Use LoopTools to compute the loop functions]'
        'mass-error-check[Check mass eigenvalue precision]'
        'meta[Create model classes]'
        'shared-libs[Create shared libraries]'
        'silent[Suppress all command line output]'
        'sqlite[Enable SQLite]'
        'static[Create statically linked executables]'
        'threads[Enable multi-threading]'
        'tsil[Enable TSIL]'
        'verbose[Verbose messages]'
        'feynarts[Enable usage of FeynArts]'
        'formcalc[Enable usage of FormCalc]'
    )

    for feature in $features; do
        args+=(
            "(--disable-${feature%%\[*})--enable-${feature}"
            "(--enable-${feature%%\[*})--disable-${feature}"
        )
    done

    _arguments -S \
        $args \
        '--with-addons=-[Comma separated list of addons to be build]:addons:_flexiblesusy_addons' \
        '--with-boost-libdir=-[Path to search for BOOST libraries]:directory:_flexiblesusy_files -/' \
        '--with-boost-incdir=-[Path to search for BOOST headers]:directory:_flexiblesusy_files -/' \
        '--with-collier-libdir=-[Path to search for COLLIER libraries]:directory:_flexiblesusy_files -/' \
        '--with-collier-incdir=-[Path to search for COLLIER modules]:directory:_flexiblesusy_files -/' \
        '--with-cppflags=-[C preprocessor flags]:flags: ' \
        '--with-cxx=-[C++ compiler]:compiler:_command_names -e' \
        '--with-cxxflags=-[C++ compiler flags]:flags: ' \
        '--with-eigen-incdir=-[Path to search for Eigen headers]:directory:_flexiblesusy_files -/' \
        '--with-fc=-[Fortran compiler]:compiler:_command_names -e' \
        '--with-fflags=-[Fortran compiler flags]:flags: ' \
        '--with-flibs=-[Fortran libraries]:libraries: ' \
        '--with-gm2calc-libdir=-[Path to search for GM2Calc library]:directory:_flexiblesusy_files -/' \
        '--with-gm2calc-incdir=-[Path to search for GM2Calc header]:directory:_flexiblesusy_files -/' \
        '--with-gsl-config=-[Path to gsl-config]:command:_command_names -e' \
        '--with-higgstools-libdir=-[Path to search for HiggsTools library]:directory:_flexiblesusy_files -/' \
        '--with-higgstools-incdir=-[Path to search for HiggsTools header]:directory:_flexiblesusy_files -/' \
        '--with-himalaya-libdir=-[Path to search for Himalaya library]:directory:_flexiblesusy_files -/' \
        '--with-himalaya-incdir=-[Path to search for Himalaya header]:directory:_flexiblesusy_files -/' \
        '--with-install-dir=-[Path to source code installation directory]:directory:_flexiblesusy_files -/' \
        '--with-lilith=-[Path to Lilith]:directory:_flexiblesusy_files -/' \
        '--with-llflags=-[Preprocessor flags for compiling with LibraryLink]:flags: ' \
        '--with-lllibs=-[Linker flags for linking LibraryLink libraries]:flags: ' \
        '--with-loop-libraries=-[Comma separated list of loop libraries to be enabled]:loop libraries:_values -s , "loop library" collier looptools fflite' \
        '--with-looptools-libdir=-[Path to search for LoopTools libraries]:directory:_flexiblesusy_files -/' \
        '--with-looptools-incdir=-[Path to search for LoopTools headers]:directory:_flexiblesusy_files -/' \
        '--with-math-cmd=-[Mathematica kernel]:command:_command_names -e' \
        '--with-models=-[Comma separated list of models to be build]:models:_flexiblesusy_created_models' \
        '--with-optional-modules=-[Load optional FlexibleSUSY modules]:modules:_values -s , module examples test' \
        '--with-pthread-libdir=-[Path to search for pthread library]:directory:_flexiblesusy_files -/' \
        '--with-python=-[Python interpreter to build the Python modules of the models for]:command:_command_names -e' \
        '--with-shared-ldflags=-[Linker flags for shared library linking]:flags: ' \
        '--with-shared-ldlibs=-[Additional libraries to be linked during shared library linking]:libraries: ' \
        '--with-shared-lib-cmd=-[Command to make shared library]:command:_command_names -e' \
        '--with-shared-lib-ext=-[File extension for shared library]:extension: ' \
        '--with-solvers=-[Comma separated list of RG solver algorithms]:solvers:_values -s , solver all two_scale lattice semi_analytic shooting' \
        '--with-sqlite-libdir=-[Path to search for SQLite library]:directory:_flexiblesusy_files -/' \
        '--with-sqlite-incdir=-[Path to search for SQLite header]:directory:_flexiblesusy_files -/' \
        '--with-static-ldflags=-[Linker flags for static linking]:flags: ' \
        '--with-static-ldlibs=-[Additional libraries to be linked during static linking]:libraries: ' \
        '--with-static-lib-cmd=-[Command to make static library]:command:_command_names -e' \
        '--with-static-lib-ext=-[File extension for static library]:extension: ' \
        '--with-tsil-libdir=-[Path to search for TSIL library]:directory:_flexiblesusy_files -/' \
        '--with-tsil-incdir=-[Path to search for TSIL header]:directory:_flexiblesusy_files -/' \
        '(- *)--help[Print the help message]' \
        '(- *)--version[Print version information]'
}

_flexiblesusy_createmodel()
{
    _arguments -S \
        '--name=-[Name of the FlexibleSUSY model]:model name: ' \
        '--model-file=-[Name of the model whose model file should be used]:model:_flexiblesusy_available_models' \
        '--sarah-model=-[Name of the SARAH model]:SARAH model: ' \
        '--with-math-cmd=-[Mathematica kernel]:command:_command_names -e' \
        '(--disable-debug)--enable-debug[Debug messages from the meta code]' \
        '(--enable-debug)--disable-debug[No debug messages from the meta code]' \
        '--force[Overwrite existing model files]' \
        '(- *)--help[Print the help message]'
}

_flexiblesusy_createaddon()
{
    _arguments -S \
        '--name=-[Name of the FlexibleSUSY addon]:addon name: ' \
        '--force[Overwrite existing addon files]' \
        '(- *)--help[Print the help message]'
}

_flexiblesusy_run_spectrum_generator()
{
    _arguments -S \
        '--slha-input-file=-[SLHA input file]:file:_flexiblesusy_files' \
        '--slha-output-file=-[SLHA output file]:file:_flexiblesusy_files' \
        '--spectrum-output-file=-[File to write spectrum to]:file:_flexiblesusy_files' \
        '--database-output-file=-[SQLite database file to write parameter point to]:file:_flexiblesusy_files' \
        '--higgsbounds-dataset=-[Location of HiggsBounds dataset]:directory:_flexiblesusy_files -/' \
        '--higgssignals-dataset=-[Location of HiggsSignals dataset]:directory:_flexiblesusy_files -/' \
        '--lilith-db=-[Lilith database]:file:_flexiblesusy_files' \
        '--rgflow-output-file=-[File to write rgflow to]:file:_flexiblesusy_files' \
        '--build-info[Print build information]' \
        '--model-info[Print model information]' \
        '(- *)--help[Print the help message]' \
        '(- *)--version[Print program version]'
}

echo "installing zsh completion for configure"
compdef _flexiblesusy_configure configure ./configure

echo "installing zsh completion for createmodel"
compdef _flexiblesusy_createmodel createmodel ./createmodel

echo "installing zsh completion for createaddon"
compdef _flexiblesusy_createaddon createaddon ./createaddon

echo "installing zsh completion for spectrum generators"
compdef -p _flexiblesusy_run_spectrum_generator '(|*/)run_*.x'
