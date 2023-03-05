_RWS_concat_or_star := rec(
           isFSA := true,
        alphabet := rec(
                type := "identifiers",
                size := 6,
              format := "dense",
               names := [a,A,b,B,t0,T0]
               ),
          states := rec(
                type := "simple",
                size := 5
               ),
           flags := ["DFA","minimized","BFS","accessible","trim"],
         initial := [1],
       accepting := [1..5],
           table := rec(
              format := "dense deterministic",
      numTransitions := 16,
         transitions := [[0,0,0,0,2,3],
                         [0,0,4,1,2,3],
                         [5,1,0,0,2,3],
                         [0,0,1,0,2,3],
                         [1,0,0,0,2,3] 
                        ]
               )
);
