_RWS_concat_or_star := rec(
           isFSA := true,
        alphabet := rec(
                type := "identifiers",
                size := 8,
              format := "dense",
               names := [a,A,b,B,c,C,t0,T0]
               ),
          states := rec(
                type := "simple",
                size := 15
               ),
           flags := ["DFA","minimized","BFS","accessible","trim"],
         initial := [1],
       accepting := [1..15],
           table := rec(
              format := "dense deterministic",
      numTransitions := 50,
         transitions := [[0,0,0,0,0,0,2,3],
                         [0,0,4,0,5,0,2,3],
                         [6,0,0,0,5,0,2,3],
                         [7,0,0,0,5,0,2,3],
                         [8,0,9,0,0,0,2,3],
                         [0,0,7,0,5,0,2,3],
                         [0,0,0,0,5,0,2,3],
                         [0,0,10,0,0,0,2,3],
                         [11,0,0,0,0,0,2,3],
                         [7,0,0,0,12,0,2,3],
                         [0,0,0,0,13,0,2,3],
                         [14,0,0,0,0,0,2,3],
                         [0,0,9,0,0,0,2,3],
                         [0,0,15,0,0,0,2,3],
                         [0,0,0,0,12,0,2,3] 
                        ]
               )
);
