_RWS_concat_or_or_star := rec(
           isFSA := true,
        alphabet := rec(
                type := "identifiers",
                size := 12,
              format := "dense",
               names := [a,A,b,B,c,C,t0,T0,t1,T1,t2,T2]
               ),
          states := rec(
                type := "simple",
                size := 13
               ),
           flags := ["DFA","minimized","BFS","accessible","trim"],
         initial := [1],
       accepting := [1..13],
           table := rec(
              format := "dense deterministic",
      numTransitions := 105,
         transitions := [[0,0,0,0,0,0,2,2,3,3,4,4],
                         [0,0,5,6,7,1,2,2,3,3,4,4],
                         [8,0,0,0,7,1,2,2,3,3,4,4],
                         [9,0,10,11,0,0,2,2,3,3,4,4],
                         [0,0,12,0,7,1,2,2,3,3,4,4],
                         [0,0,0,8,7,1,2,2,3,3,4,4],
                         [0,0,0,0,1,0,2,2,3,3,4,4],
                         [0,0,0,0,7,1,2,2,3,3,4,4],
                         [0,0,10,11,0,0,2,2,3,3,4,4],
                         [0,0,13,0,0,0,2,2,3,3,4,4],
                         [0,0,0,1,0,0,2,2,3,3,4,4],
                         [0,0,8,0,7,1,2,2,3,3,4,4],
                         [0,0,1,0,0,0,2,2,3,3,4,4] 
                        ]
               )
);
