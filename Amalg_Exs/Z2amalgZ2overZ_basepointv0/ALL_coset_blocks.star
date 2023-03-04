_RWS_concat_or_star := rec(
           isFSA := true,
        alphabet := rec(
                type := "identifiers",
                size := 10,
              format := "dense",
               names := [a,A,b,B,c,C,d,D,t0,T0]
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
         transitions := [[0,0,0,0,0,0,0,0,2,3],
                         [0,0,0,0,0,0,4,5,2,3],
                         [0,0,6,7,0,0,0,0,2,3],
                         [0,0,0,0,0,0,4,0,2,3],
                         [0,0,0,0,0,0,0,5,2,3],
                         [0,0,6,0,0,0,0,0,2,3],
                         [0,0,0,7,0,0,0,0,2,3] 
                        ]
               )
);
