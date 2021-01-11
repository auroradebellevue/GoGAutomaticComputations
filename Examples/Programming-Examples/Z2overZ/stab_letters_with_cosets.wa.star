_RWS_concat_or_star := rec(
           isFSA := true,
        alphabet := rec(
                type := "identifiers",
                size := 6,
              format := "dense",
               names := [b,B,a,A,t,T]
               ),
          states := rec(
                type := "simple",
                size := 7
               ),
           flags := ["DFA","minimized","BFS","accessible","trim"],
         initial := [1],
       accepting := [1..7],
           table := rec(
              format := "dense deterministic",
      numTransitions := 22,
         transitions := [[0,0,0,0,2,3],
                         [4,5,0,0,2,3],
                         [0,0,6,7,2,3],
                         [4,0,0,0,2,3],
                         [0,5,0,0,2,3],
                         [0,0,6,0,2,3],
                         [0,0,0,7,2,3] 
                        ]
               )
);
