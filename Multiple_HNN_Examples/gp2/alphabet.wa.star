_RWS.wa_star := rec(
           isFSA := true,
        alphabet := rec(
                type := "identifiers",
                size := 12,
              format := "dense",
               names := [a,A,b,B,c,C,t0,T0,t1,T1,t2,T2]
               ),
          states := rec(
                type := "simple",
                size := 1
               ),
           flags := ["DFA","minimized","BFS","accessible","trim"],
         initial := [1],
       accepting := [1],
           table := rec(
              format := "dense deterministic",
      numTransitions := 12,
         transitions := [[1,1,1,1,1,1,1,1,1,1,1,1] 
                        ]
               )
);
