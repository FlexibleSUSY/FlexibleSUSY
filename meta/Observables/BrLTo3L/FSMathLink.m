Begin@"FSMathLink`Private`";

PutObservable[obs:FlexibleSUSYObservable`BrLTo3L[_, _ -> {_, _, _}, _, _], type_, link_String, heads_:{}] :=
    PutFirstEntryOfObservable[obs, link];

End[];
