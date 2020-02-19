==================
Add new observable
==================

FlexibleSUSY provides a set of observables to be calculated in
addition to the model spectrum itself.

At some point, to make deeper analysis of a model one would like to
add some new, or modify existing observables. This file can be used then
as a reference to see generic stucture of this implementation.

Mathematica part
````````````````
Open file ``meta/Observables.m`` and do there the following steps:

1) Add a ``Mathematica`` name for new observable to ``FSObservables``.

2) Add new definition for a ``GetObservableName`` and ``GetObservableDescription`` functions.

 Note: ``GetObservableName`` specifies a c++ name for a function (further
 ``<name_of_observable>`` is used), so it should give a unique valid one.
 There is a convention, that these name on c++ level
 should be underline separated, with usage of small letters.

 ``GetObservableDescription`` should be a ``string`` with
 desription of some useful for user information.

3) Modify ``GetRequestedObservables`` by imposing some restrictions for
   a new observable.

 Note: one possibly would like to move **all** observable restrictions here.

4) Add ``GetObservableType`` for a new observable.

 Note: available types are given in ``meta/CConversion.m`` file.

5) Add ``CalculateObservable``.

 Note: it should return a valid c++ command. On ``Mathematica`` level it
 is usually written in the following way::

   structName <>
   ".<NameOfObservable><0|1>(<arguments>) = " <>
   "<ActualCallOfObservable>;"

 ``<NameOfObservable><0,1>(<arguments>)`` is a c++ preprocessor macro
 which is converted to the name of desired observable ``<name_of_observable>``.

 * ``<NameOfObservable>`` is a beginning of the this macro (should be unique,
   see ``templates/observables.cpp`` for actual set of used definitions).
 * ``<0,1>`` is either ``0`` or ``1``. Internally there is a convention, that
   observable should accept fields with and without generation indices.
   ``<0,1>`` is replaced by ``1`` for the first case and by ``0`` for the second one.
 * ``<arguments>`` is a comma separated set of macros arguments which specify the
   name of a observable.

 There is an interesting question: why is it done in this way? Because actually,
 this macros should return **exactly** the same expression as the one by
 ``GetObservableName``. Why not just call this Mathematica function?

 ``<ActualCallOfObservable>`` is a function call which actually calculates an
 observable. Usually, it is called in the following way::

  <namespace_of_this_observable>::<calculate_function>(<args>)

 where ``<args>`` can differ from ``<arguments>`` of c++ preprocessor macro call.

Open file ``meta/FlexibleSUSY.m``

1) Inside ``BeginPackage`` add ``"<NewMetaFile>`"``. This corresponds to a new
   file to be created later ``meta/<NewMetaFile>.m``.

2) Create somewhere (close to analogous definitions) a new function
   ``Write<NewObservable>Class>``, which should (usual input give in the next item):

* convert (modified somehow)
  ``Observables`GetRequestedObservables[extraSLHAOutputBlocks]``
  into some c++ code;

* read a content of files (to be created later in this
  manual) ``template/<new_observable>.{hpp.in,cpp.in}``;

* replace some tokens;

* save the result to ``models/<ModelName>/<new_observable>.{hpp,cpp}``;

* return ``List`` of ``List`` of field names, which represent required by
  calculation vertices.


3) Inside ``MakeFlexibleSUSY`` function add ``<returnedVertices>`` to a list of
   local symbols and::

     (* OBSERVABLE: <some description> *)
     Print["Creating <> class ..."];
     <returnedVertices> =
        Write<NewObservable>Class[
           <main input>,
           {
              {
                 FileNameJoin@{$flexiblesusyTemplateDir, #<>".hpp.in"},
                 FileNameJoin@{FSOutputDir, FlexibleSUSY`FSModelName<>"_"<>#<>".hpp"}
              },
              {
                 FileNameJoin@{$flexiblesusyTemplateDir, #<>".cpp.in"},
                 FileNameJoin@{FSOutputDir, FlexibleSUSY`FSModelName<>"_"<>#<>".cpp"}
              }
           } & ["<new_observable>"]
     ];

   to the place where other classes are created.

   A little bit further add ``<returnedVertices>`` to a ``Join`` of first argument of
   ``WriteCXXDiagramClass``.

Create ``meta/<NewMetaFile>.m``

1) It should make some model specific calculation of an observable.

Now FlexibleSUSY knows on meta level how to write ``observables.{hpp,cpp}`` files.
This is not enough, though. One also need to specify explicitly, how new
observable should be formatted in ``SHLA`` file.

To do this, open ``meta/WriteOut.m`` file and edit ``WriteSLHABlockEntry``
function by adding information about new variable.

 Note: there are several overload ``WriteSLHABlockEntry`` functions. In the one
 which checks for ``?IsObservable`` one adds generic behavior of this function
 to a specific for an observable.

Type of observable is not supported
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In file ``src/shla_io.hpp`` there are several ``#define FORMAT_<SOMETHING>``
macros. One can extend them to include some new formatting for not yet
included types. Or, what is a little bit better, one can use them as building
blocks to construct some additional, more complex and flexible output.



C++ part
````````

Open file ``templates/observables.cpp``

1) Add include preprocessor command for header of your observable::

     #include "@ModelName@_<new_observable>.hpp"

2) Add define preprocessor command in the form::

     #define  <NameOfObservable><0,1>(<arguments>) <name_of_observable>

   which will generate the name of observable.

3) Create files ``template/<new_observable>.{hpp.in,cpp.in}``, where actual c++
   calculations are done. The only requirements is that function
   ``flexiblessusy::<namespace_of_this_observable>::<calculate_function>(<args>)`` exists.

Automatization part
```````````````````
Open file ``templates/module.mk`` and add::

   $(DIR)/<new_observable>.hpp.in \
   $(DIR)/<new_observable>.cpp.in \

to a ``BASE_TEMPLATES`` variable.

Open file ``templates/module.mk.in`` and add::

   $(DIR)/@CLASSNAME@_<new_observable>.cpp \

to a ``LIB@CLASSNAME@_SRC`` variable. Also add::

   $(DIR)/@CLASSNAME@_<new_observable>.hpp \

to a ``LIB@CLASSNAME@_HDR`` variable.

Open file ``meta/module.mk`` and add::

   $(DIR)/<NewMetaFile>.m \

to ``META_SRC`` variable.

 Note: this allows correct code building if ``.m`` file was changed.

Documentation part
``````````````````
