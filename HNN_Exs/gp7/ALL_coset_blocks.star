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
                size := 19
               ),
           flags := ["DFA","minimized","BFS","accessible","trim"],
         initial := [1],
       accepting := [1..19],
           table := rec(
              format := "dense deterministic",
      numTransitions := 74,
         transitions := [[0,0,0,0,2,3],
                         [0,0,4,5,2,3],
                         [6,7,0,0,2,3],
                         [8,9,4,0,2,3],
                         [10,11,0,5,2,3],
                         [6,0,12,13,2,3],
                         [0,7,14,15,2,3],
                         [16,0,0,0,2,3],
                         [0,9,4,0,2,3],
                         [10,0,0,5,2,3],
                         [0,17,0,0,2,3],
                         [0,0,18,0,2,3],
                         [6,0,0,13,2,3],
                         [0,7,14,0,2,3],
                         [0,0,0,19,2,3],
                         [16,0,12,0,2,3],
                         [0,17,0,15,2,3],
                         [8,0,18,0,2,3],
                         [0,11,0,19,2,3] 
                        ]
               )
);
