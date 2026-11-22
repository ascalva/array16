<!-- Written by the `report` step of array16.lab. -->

# The verdicts

14 of 14 shots are on file.

| Integrity |  |
| --- | --- |
| shots | 14 |
| measured | 14 |
| predictions match the code | True |
| records made before registration | 0 |
| records from other predictions | 0 |
| records from uncommitted code | 0 |
| commits | ['3e2149f6fc2c'] |
| benches | ['twin-r1'] |

| Criterion | Passed | Judged |
| --- | --- | --- |
| E2-21 well | 2 | 2 |
| E2-21 well control | 1 | 1 |
| slot depth | 1 | 1 |
| E2-21 hold | 1 | 2 |
| E2-21 senses | 3 | 3 |
| E2-21 hold control | 1 | 1 |
| E2-22 static limit | 2 | 3 |
| E2-22 price of the knot | 1 | 2 |
| slot secular frequency | 1 | 1 |
| must fail control | 1 | 1 |
| must fail | 2 | 2 |
| E2-21, every check | 8 | 9 |
| E2-22, every check | 3 | 5 |
| slot, every check | 2 | 2 |
| must fail, every check | 3 | 3 |

## Failed: E2-21 hold

- E2-21 hold, R1, realisation 0: held to 9.344 eV (with the wave 9.344 eV, against it 9.344 eV; the first loss at 9.229 eV) against the static limit 9.307 eV (band 5.677 to 8.842)

## Failed: E2-22 static limit

- E2-22 hold, flat ring 318, realisation 0: held to 32.25 eV (with the wave 32.25 eV, against it 32.25 eV; the first loss at 32.65 eV) against the static limit 32.43 eV (band 19.78 to 30.81)

## Failed: E2-22 price of the knot

- E2-22 hold, flat ring 318, realisation 0: 32.25 eV over T915-S's 9.996 eV: 3.226 against the static 2.786 (band 2.368 to 3.204)

## Failed: E2-21, every check

- E2-21 hold, R1, realisation 0: held to 9.344 eV (with the wave 9.344 eV, against it 9.344 eV; the first loss at 9.229 eV) against the static limit 9.307 eV (band 5.677 to 8.842)

## Failed: E2-22, every check

- E2-22 hold, flat ring 318, realisation 0: held to 32.25 eV (with the wave 32.25 eV, against it 32.25 eV; the first loss at 32.65 eV) against the static limit 32.43 eV (band 19.78 to 30.81)
- E2-22 hold, flat ring 318, realisation 0: 32.25 eV over T915-S's 9.996 eV: 3.226 against the static 2.786 (band 2.368 to 3.204)

## Every comparison

| Criterion | Shot | Result | Numbers |
| --- | --- | --- | --- |
| E2-21 well | E2-21 well, R1, realisation 0 | pass | first lost at 0.01179 eV against the axis model's 0.01179 eV (band 0.01108 to 0.01254); the robust threshold 0.03367 eV |
| E2-21 well | E2-21 well, R1, on the bench | pass | first lost at 0.01547 eV against the axis model's 0.01179 eV (band 0.01081 to 0.02268); the robust threshold 0.01728 eV |
| E2-21 well control | E2-21 well, R1, realisation 0 at q 0.3 | pass | first lost at 0.03138 eV against the axis model's 0.03138 eV (band 0.0295 to 0.03338); the robust threshold 0.03256 eV |
| slot depth | slot well, straight slot guide, realisation 0 | pass | escape at 0.2074 eV against the depth 0.2294 eV (band 0.1835 to 0.2867) |
| E2-21 hold | E2-21 hold, R1, realisation 0 | **FAIL** | held to 9.344 eV (with the wave 9.344 eV, against it 9.344 eV; the first loss at 9.229 eV) against the static limit 9.307 eV (band 5.677 to 8.842) |
| E2-21 senses | E2-21 hold, R1, realisation 0 | pass | the senses' limits 0.0 scan steps apart, 2 allowed |
| E2-21 hold | E2-21 hold, R1, on the bench | pass | held to 8.465 eV (with the wave 7.86 eV, against it 9.116 eV; the first loss at 7.764 eV) against the static limit 9.307 eV (band 4.69 to 9.3) |
| E2-21 senses | E2-21 hold, R1, on the bench | pass | |sqrt W- - sqrt W+| = 0.216 sqrt(eV) against Earth's most and two steps, 0.670 |
| E2-21 hold control | E2-21 hold, R1, realisation 0 at q 0.3 | pass | held to 3.269 eV (with the wave 3.309 eV, against it 3.228 eV; the first loss at 3.269 eV) against the static limit 3.28 eV (band 2.788 to 3.444) |
| E2-21 senses | E2-21 hold, R1, realisation 0 at q 0.3 | pass | the senses' limits 1.0 scan steps apart, 2 allowed |
| E2-22 static limit | E2-22 hold, T915-S, realisation 0 | pass | held to 9.996 eV (with the wave 9.282 eV, against it 10.76 eV; the first loss at 9.169 eV) against the static limit 11.64 eV (band 7.102 to 11.06) |
| E2-22 static limit | E2-22 hold, flat ring 318, realisation 0 | **FAIL** | held to 32.25 eV (with the wave 32.25 eV, against it 32.25 eV; the first loss at 32.65 eV) against the static limit 32.43 eV (band 19.78 to 30.81) |
| E2-22 static limit | E2-22 hold, helix (1, 3), realisation 0 | pass | held to 16.1 eV (with the wave 16.1 eV, against it 16.1 eV; the first loss at 16.3 eV) against the static limit 18.45 eV (band 11.26 to 17.53) |
| E2-22 price of the knot | E2-22 hold, flat ring 318, realisation 0 | **FAIL** | 32.25 eV over T915-S's 9.996 eV: 3.226 against the static 2.786 (band 2.368 to 3.204) |
| E2-22 price of the knot | E2-22 hold, helix (1, 3), realisation 0 | pass | 16.1 eV over T915-S's 9.996 eV: 1.611 against the static 1.585 (band 1.347 to 1.822) |
| slot secular frequency | slot beta, straight slot guide, realisation 0 | pass | beta 0.326070 to 0.326078 against Mathieu's 0.326034 (at most 1.36e-04 off, 0.002 allowed); 4 of 4 held |
| must fail control | must fail control, R1, on the bench | pass | q 0.51 at 50.0 V: 32 of 32 held, 0 lost to the guide |
| must fail | must fail, R1, on the bench at a 1 eV well | pass | q 1.66 at 162.4 V (theory 1.64): 0 of 32 held, the last lost after 2.2 RF periods; the control held |
| must fail | must fail, R1, on the bench at q 1.89 | pass | q 1.91 at 187.0 V (theory 1.89): 0 of 32 held, the last lost after 1.9 RF periods; the control held |
