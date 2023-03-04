_RWS.wa_star := rec(
           isFSA := true,
        alphabet := rec(
                type := "identifiers",
                size := 6,
              format := "dense",
               names := [a,A,b,B,t0,T0]
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
      numTransitions := 6,
         transitions := [[1,1,1,1,1,1] 
                        ]
               )
);
