==================
Add new observable
==================

FlexibleSUSY provides a set of observables to be calculated in
addition to the model spectrum itself.

At some point, to make deeper a analysis of a model one would like to
add some new, or modify existing observables. This file can be used
as a reference to see a generic stucture of this implementation.

Mathematica part
````````````````

We need to specify some FlexibleSUSY-level C++/Mathematica relations for
a new observable. It can be done with the help of ``NPointFunctions`` module as
follows.

There is a file ``meta/Observables.m`` file, which contains some generic
information about each observable. For the sake of code orthogonality the parts,
responsible for different observables, can be moved into different files.
The convention for this separation:

  Information is stored in
  ``meta/NPointFunctions/<NameOfObservable>/observable.m`` files. And it is
  loaded only if ``FeynArts`` and ``FormCalc`` are allowed to be used at
  ``configure`` stage. ``<NameOfObservable>`` exactly corresponds to
  the Mathematica name of an observable and its syntax is represented in by
  the name (every new word with a large letter - usual Mathematica convention).

So the steps are the following:

Create a ``meta/NPointFunctions/<NameOfObservable>/type.m`` file.

* **Hint:** one can copy some form existing ones in order to simplify
  the modification.

1) This file should contain some "global" static things, such as C++ convention
   for a namespace of an observable, arguments allowed, etc.
   See the content of existing files to get inspiration and pay attention to
   some syntactic sugar of context and package definitions.

Create a ``meta/NPointFunctions/<NameOfObservable>/observable.m`` file.

0) Link it to a corresponding ``type.m`` file, such that it is loaded
   automatically.

1) Add new definition for a ``GetObservableName`` and ``GetObservableDescription``
   functions.

   * ``GetObservableName`` specifies a C++ name for a function (further
     ``<name_of_observable>`` is used), so it should give a *unique* valid one.
     **Hint:** One should use here *explicit* generation numbers.
     There is a convention, that these names on C++ level should be underline
     separated, with usage of small letters (like in the pattern above).

   * ``GetObservableDescription`` should be a ``string`` with desription of
     some useful information, related to observable.

2) Add ``GetObservableType`` for a new observable.

   * **Hint:** Available types are given in ``meta/CConversion.m`` file.

3) Add ``CalculateObservable``.

   * **Hint:** It should return a valid C++ command. On ``Mathematica`` level it
     is usually written in the following way::

         structName <> ".<name_of_observable> = <ActualCallOfObservable>;"

     where ``<ActualCallOfObservable>`` is a function call which actually calculates an
     observable. Usually, it is called in the following way::

         <model>_<observable_namespace>::<calculate_function>(<args>)

Create a ``meta/NPointFunctions/<NameOfObservable>/write.m`` file. It will
specify the SLHA/FLHA outputs of and observable and will be loaded by ``WriteOut``
(if ``FeynArts`` and ``FormCalc`` are avaliable).

1) It should have the function ``<NameOfObservable>`write`` function with the
   following arguments::

      _String,  (* name of Les Houches block *)
      _,        (* observable with arguments *)
      _Integer, (* number of block entry *)
      _String   (* comment for block entry *)

   which generates relevant C++ code, which will specify the output (**hint:**
   see existing files for inspiration).

Create a ``meta/NPointFunctions/<NameOfObservable>/class.m`` file with
``Write<NewObservable>Class`` function, which should:

1) Have arguments of the following format (with arbitrary names)
   ``_List, {{_?FileExistsQ, _String}..}``.

2) Read a content of files (to be created some steps later)
   ``template/<observable_namespace>.{hpp.in,cpp.in}``, modify them and save
   into ``models/<model>/<observable_namespace>.{hpp,cpp}``.

3) Return two ``List`` with some fields and some vertices.

**Optional**: Open file ``meta/FlexibleSUSY.m`` and modify ``MakeFlexibleSUSY``
body in case if some additional actions are required. For example, first
returned values of ``Write<NewObservable>Class`` are saved in ``FieldsNPF@"<NewObservable>"``.
The second ones - in ``VerticesNPF@"<NewObservable>"``. They can be used to
make some additional actions.

Type of observable is not supported
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In file ``src/shla_io.hpp`` there are several ``#define FORMAT_<SOMETHING>``
macros. One can extend them to include some new formatting for not yet
included types.


C++ part
````````

Open file ``templates/observables.cpp.in``

1) Add include preprocessor command for header of your observable::

     #include "@ModelName@_<observable_namespace>.hpp"

2) Create files ``template/<observable_namespace>.{hpp.in,cpp.in}``, where actual C++
   calculations are done. The only requirement is that function
   ``flexiblessusy::<observable_namespace>::<calculate_function>(<args>)`` exists.

Automatization part
```````````````````
Open file ``templates/module.mk`` and add::

   $(DIR)/<observable_namespace>.hpp.in \
   $(DIR)/<observable_namespace>.cpp.in \

to a ``BASE_TEMPLATES`` variable.

Open file ``templates/module.mk.in`` and add::

   $(DIR)/@CLASSNAME@_<new_observable>.cpp \

to a ``LIB@CLASSNAME@_SRC`` variable. Also add::

   $(DIR)/@CLASSNAME@_<new_observable>.hpp \

to a ``LIB@CLASSNAME@_HDR`` variable.

Documentation part
``````````````````
