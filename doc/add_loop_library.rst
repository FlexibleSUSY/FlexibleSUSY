===============
Add new library
===============

FlexibleSUSY provides a loop library interface, which allows to setup several
loop libraries once and then choose any of them at run time.

If one has a wish to add another library to a set of existing ones then this
document can be used in order to simplify required operations.

Main steps and definitions
``````````````````````````

All procedure can be divided into the following set:

* Create C++ code of the library
* Make ``./configure`` automatically work with new library
* Modify documentation

The following abbreviations are used:

``<Newlibrary>`` - name of loop library for a ``enum class`` inside
``loop_library.hpp`` and a class name for a library with the following hierarchy
``flexiblesusy::looplibrary::<Newlibrary>``.
First letter is capital, all others are small.

``<NEWLIBRARY>`` - name of loop library used for C++ preprocessor commands and
``configure`` script. All letters are capital.

``<newlibrary>`` - name of loop library used for user input for ``configure``
script and for part of C++ code names.
All letters are small.

``<N>`` - ``int`` number for a new loop library. Should be unique.

Create C++ code of the library
``````````````````````````````

We need let ``FlexibleSUSY`` to know that there is a new library. Information
about this is stored inside ``src/loop_libraries/loop_library.{hpp,cpp}`` files.


* Open a ``.hpp`` file and add a name of a new library into ``enum class Library``.
  Currently there is a convention that the name has a form ``<Newlibrary>``.

* Open a ``.cpp`` file and add in the head of the file after similar expressions
  the following::

    #ifdef ENABLE_<NEWLIBRARY>
    #include "library_<newlibrary>.hpp"
    #define <NEWLIBRARY>_INFO ", <N> (=<newlibrary>)"
    #else
    #define <NEWLIBRARY>_INFO
    #endif // ENABLE_<NEWLIBRARY>

  Add a new entry to a switch statement right before ``default:``::

    #ifdef ENABLE_<NEWLIBRARY>
    case <N> : Loop_library::lib_ = std::make_unique<looplibrary::<Newlibrary>>();
             Loop_library::type_ = Loop_library::Library::<Newlibrary>;
             break;
    #endif // ENABLE_<NEWLIBRARY>

  Add ``<NEWLIBRARY>_INFO`` to ``invalid_argument()`` function.

Note that ``ENABLE_<NEWLIBRARY>`` is used here. Whether it is enabled or disabled
is (actually will be) known from ``config/config.h`` file, which (actually
``config/config.h.in`` file) will be modified
later in this manual to include this information.

* Create ``src/loop_libraries/library_<newlibrary>.{hpp,cpp}`` files. One can
  look at already existing libraries in directory ``src/loop_libraries`` to get
  knowledge about used names. For example, ``.hpp`` file require extremely few
  modification thanks to ``boost/preprocessor`` macros.

Make configure routines
```````````````````````

In the previous step we created a main C++ code for a new library. Now we need
let ``FlexibleSUSY`` to know how to be compiled with this library.

* In file ``src/module.mk`` add after analogous definitions for other loop libraries
  (anchor ``# loop library #`` can be used as a starting point) the following
  lines::

    ifeq ($(ENABLE_<NEWLIBRARY>),yes)
    LIBFLEXI_SRC += \
        $(DIR)/loop_libraries/library_<newlibrary>.cpp
    LIBFLEXI_HDR += \
        $(DIR)/loop_libraries/library_<newlibrary>.hpp
    endif

  One can look here at the implementation of additional commands for COLLIER library
  to get an idea of how to add some additional dependencies here.

We added check for a new yet undefined variable ``ENABLE_<NEWLIBRARY>``. Now we need
to modify ``configure`` file mainly in order to setup correctly this variable.

Open ``configure`` file and there do the following steps:

* Add line::

    enable_<newlibrary>="no"

  in the beginning of the documents right after analogous variables.

* After ``# corresponding preprocessor define statements`` add::

    DEFINE_ENABLE_<NEWLIBRARY>="#undef ENABLE_<NEWLIBRARY>"

  and after these DEFINE statements add::

    <newlibrary>_lib_dir=""
    <newlibrary>_inc_dir=""

* Inside function ``enable_defines()`` add::

    if [ "x$enable_<newlibrary>" = "xyes" ]; then
       DEFINE_ENABLE_<NEWLIBRARY>="#define ENABLE_<NEWLIBRARY> 1"
       message "Enabling usage of <newlibrary>"
       logmsg "   ${DEFINE_ENABLE_<NEWLIBRARY>}"
    else
       enable_<newlibrary>="no"
       DEFINE_ENABLE_<NEWLIBRARY>="#undef ENABLE_<NEWLIBRARY>"
       logmsg "Disabling usage of <newlibrary>"
       logmsg "   ${DEFINE_ENABLE_<NEWLIBRARY>}"
    fi
    logmsg "   ${DEFINE_ENABLE_<NEWLIBRARY>}"

  This variables will go to ``config/config.h.in`` afterwards.

* Add inside ``add_metaflags()``::

    test $enable_<newlibrary> = 'yes' && lib_="${lib_}, FS<NEWLIBRARY>"

* Add inside ``replace_markers()``::

    -e "s|@ENABLE_<NEWLIBRARY>@|$enable_<newlibrary>|" \
    -e "s|@<NEWLIBRARY>LIBS@|$<NEWLIBRARY>LIBS|"       \
    -e "s|@<NEWLIBRARY>FLAGS@|$<NEWLIBRARY>FLAGS|"     \

  This will go to (inside ``config/`` subdirectory) ``flexiblesusy-config.in``
  (``Makefile.standalone.in``, ``Makefile.customized-betas.in``,
  ``Makefile.tower.in`` go to null) and to
  ``Makefile.in``, which is the most important one.

* Add before ``< $CONFIGHDR_TMPL > $CONFIGHDR``::

    -e "s|@DEFINE_ENABLE_<NEWLIBRARY>@|$DEFINE_ENABLE_<NEWLIBRARY>|" \
    -e "s|@<NEWLIBRARY>FLAGS@|$<NEWLIBRARY>FLAGS|"                   \
    -e "s|@<NEWLIBRARY>LIBS@|$<NEWLIBRARY>LIBS|"                     \

* Add to the part of the code which checks cmd arguments (``if test $# -gt 0 ; then``)::

    --with-<newlibrary>-libdir=*) <newlibrary>_lib_dir=$optarg ;;
    --with-<newlibrary>-incdir=*) <newlibrary>_inc_dir=$optarg ;;

* Add ``check_<newlibrary>()``, ``check_<newlibrary>_incl()``,
  ``check_<newlibrary>_libs()``
  functions with desired behavior and structure similar to already existing ones and
  add them to function evaluation sequence (large set of functions called one after
  other).

  Note that first one usually checks settings in ``.pc`` files and runs other two.
  ``_incl()`` checks include directories and defines ``<NEWLIBRARY>FLAGS`` variable.
  ``_libs()`` checks library directories and defines ``<NEWLIBRARY>LIBS`` variable.

* Modify ``check_looplibrary()`` in the way similar to existing one. Note that
  currently old options are checked first and then there is a check for consistent
  choice of libraries. Then ``check_<library>`` functions run and modification of
  ``LOOPFUNCLIBS`` and ``LOOPFUNCFLAGS`` is performed.

Now goes a chain of changes in some additional files which will influence compilation
itself more directly.

* Open file ``config/Makefile.in`` and add after ``# Makefile`` switches::

    ENABLE_<NEWLIBRARY> := @ENABLE_<NEWLIBRARY>@

  add after ``# Variables for compilation``::

    <NEWLIBRARY>FLAGS := @<NEWLIBRARY>FLAGS@
    <NEWLIBRARY>LIBS  := @<NEWLIBRARY>LIBS@

* Open file ``config/config.h.in`` and add after ``/* Build variables */``::

    #define <NEWLIBRARY>FLAGS   "@<NEWLIBRARY>FLAGS@"
    #define <NEWLIBRARY>LIBS    "@<NEWLIBRARY>LIBS@"

  add after ``/* Switches */``::

    /* Enable <newlibrary> */
    @DEFINE_ENABLE_<NEWLIBRARY>@

* Open file ``config/Makefile.standalone.in`` and add after ``# Switches``::

    ENABLE_<NEWLIBRARY> := @ENABLE_<NEWLIBRARY>@

  add after ``showbuild``::

    @echo "ENABLE_<NEWLIBRARY> = $(ENABLE_<NEWLIBRARY>)"

  add to places where ENABLE of other libraries present::

    $(ENABLE_<NEWLIBRARY>)

* For file ``config/Makefile.tower.in`` repeat instructions for ``Makefile.standalone.in``.
* For file ``config/Makefile.customized-betas.in`` repeat instructions for ``Makefile.standalone.in``.

Note: file ``config/flexiblesusy-config.in`` could be be but was not modified by the
author of this manual.

* Open fie ``meta/FlexibleSUSY.m`` and add after ``FSLoopLibrary::usage``::

    FS<NEWLIBRARY>;

Modify documentation
````````````````````

* Inside ``configure`` change ``help()`` by adding::

   --with-<newlibrary>-libdir=    Path to search for <NEWLIBRARY> libraries
   --with-<newlibrary>-incdir=    Path to search for <NEWLIBRARY> modules

  and modifying sentence after::

   --with-loop-libraries=

* In file ``src/atom src/spectrum_generator_settings.cpp`` modify table
  "Resets all spectrum generator settings to their defaults." by adding description
  of new library

* Modify description of ``FlexibleSUSY[31]`` in ``doc/slha_input.rst`` file.

* Add description of new library inside ``README.rst`` file.
