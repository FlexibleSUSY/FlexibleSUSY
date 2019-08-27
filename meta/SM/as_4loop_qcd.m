(* QCD strong coupling threshold correction from [hep-ph/0512060].

 alphaS = strong coupling in QCD with nf flavours

 alphaSPrime = strong coupling in QCD with nl = (nf - 1) flavours

 L = Log[Q^2/mf[Q]^2] see below Eq.(27)

 Q = renormalization scale

 mf[Q] = running heavy quark mass

 *)

a4 = PolyLog[4, 1/2]

a5 = PolyLog[5, 1/2]

T543 = -8445.8046390310298

T622 = -4553.4004372195263

(* [hep-ph/0512060] Eq.(27) *)
DeltaMS4 = (
(
     134805853579559/43342154956800
    -254709337/783820800 Pi^4
    -151369/30481920 Pi^6
    -18233772727/783820800 Zeta[3]
    +151369/544320 Zeta[3]^2
    +4330717/207360 Zeta[5]
    +9869857/272160 a4
    -121/36 a5
    -2057/51840 Pi^4 Log[2]
    -9869857/6531840 Pi^2 Log[2]^2
    -121/2592 Pi^2 Log[2]^3
    +9869857/6531840 Log[2]^4
    +121/4320 Log[2]^5
    +82037/30965760 T543
    -151369/11612160 T622
)

+ nl (
    -4770941/2239488
    -541549/14929920 Pi^4
    +3645913/995328 Zeta[3]
    +115/576 Zeta[5]
    +685/5184 a4
    -685/124416 Pi^2 Log[2]^2
    +685/124416 Log[2]^4
)

+ nl^2 (
    -271883/4478976
    +167/5184 Zeta[3]
)
)

(* [hep-ph/0512060] Eq.(23) *)
dMS1 = -1/6 L

(* [hep-ph/0512060] Eq.(24) *)
dMS2 = 11/72 - 11/24 L + 1/36 L^2

(* [hep-ph/0512060] Eq.(25) *)
dMS3 = (
    + 564731/124416
    - 82043/27648 Zeta[3]
    - 955/576 L
    + 53/576 L^2
    - 1/216 L^3
    + nl (
        - 2633/31104
        + 67/576 L
        - 1/36 L^2
    )
)

(* [hep-ph/0512060] Eq.(26) *)
dMS4 = (
    DeltaMS4
    + 7391699/746496 L
    - 2529743/165888 Zeta[3] L
    + 2177/3456  L^2
    - 1883/10368  L^3
    + 1/1296 L^4
    + nl (
        - 110341/373248  L
        + 110779/82944 Zeta[3] L
        - 1483/10368 L^2
        - 127/5184 L^3
    )
    + nl^2 (
        + 6865/186624 L
        - 77/20736 L^2
        + 1/324 L^3
    )
)

as = alphaS / Pi

(* [hep-ph/0512060] Eq.(22) *)
zetag2 = (
    1
    + as^1 dMS1
    + as^2 dMS2
    + as^3 dMS3
    + as^4 dMS4
)

(* [hep-ph/0512060] Eq.(35) *)
dprimeMS1 = 1/6 L

(* [hep-ph/0512060] Eq.(36) *)
dprimeMS2 = -11/72 + 11/24 L + 1/36 L^2

(* [hep-ph/0512060] Eq.(37) *)
dprimeMS3 = (
    - 564731/124416
    + 82043/27648 Zeta[3]
    + 2645/1728 L
    + 167/576 L^2
    + 1/216 L^3
    + nl (2633/31104 - 67/576 L + 1/36 L^2)
)

(* [hep-ph/0512060] Eq.(38) *)
dprimeMS4 = (
    + 121/1728
    - DeltaMS4
    - 11093717/746496 L
    + 3022001/165888 Zeta[3] L
    + 1837/1152 L^2
    + 2909/10368 L^3
    + 1/1296 L^4
    + nl (
       + 141937/373248 L
       - 110779/82944 Zeta[3] L
       + 277/10368 L^2
       + 271/5184 L^3
    )
    + nl^2 (
       - 6865/186624 L
       + 77/20736 L^2
       - 1/324 L^3
    )
)

asprime = alphaSPrime / Pi

(* [hep-ph/0512060] Eq.(33) *)
oneOverZetag2 = (
    1
    + asprime^1 dprimeMS1
    + asprime^2 dprimeMS2
    + asprime^3 dprimeMS3
    + asprime^4 dprimeMS4
)

(* [hep-ph/0512060] Eq.(14) *)
(* alphaSPrime = alphaS zetag2 *)

(* [hep-ph/0512060] Eq.(14), solved for alpha_s(Q) *)
alphaS = alphaSPrime oneOverZetag2
