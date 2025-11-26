# Meine dollen Wellen

## Erzwungene Welle
OK 0. #ifdef für Watchpoints
OK (nicht machen )1. point source as InterpolatedPressure BC
OK 2. richtige IndicatorSphere
OK 3. Material 4 außen (check innerExtend und innerOrigin) --> in geometry.pvd nachgucken
OK 4. outdir naming

5. minimal kleinere Auflösung
OK (Hier nochmal nachschauen wie groß genau) 6. domäne ausdehnung alle Raumrichtungen ca. 3x lambda - falls zu langsam, erstmal 2xlambda
OK 6. Die Messmethode an das neue dx anpassen. dx wird mittlerweile immer neu berechnet in Abhängigkeit von der Wellenlänge
OK 7. eine kleinere Punktquelle gestalten


--------------------------
Freie Welle
OK 8aa. Logdatei schreiben, damit der Terminalinhalt mit ausgegeben wird
OK 8a. periodische Randbedingungen der freien Welle überprüfen und zusehen, dass die Domäne so lang ist, dass die Wellen komplett abschließen
OK 8. Gitterpunkte je Wellenzahl in die Ausgabe einfügen
OK 9. Viskosität neu ansetzen
10. cp analytisch sauber in das Programm programmieren, damit du weißt, ob die Werte richtig sind!
11. Überall bei den _lat Ausgaben _LU draus machen
12. gp/lamda verändert sich momentan. LDomäne und alls bleibt exakt gleich. Soll das so sein?



Erzwungene Welle
OK 8aa. Logdatei schreiben, damit der Terminalinhalt mit ausgegeben wird
OK 8. Gitterpunkte je Wellenzahl in die Ausgabe einfügen
OK9. Viskosität neu ansetzen
9. Messpunkte nicht in Abhängigkeit von lambda. Diagonal auf die Achse zusätzlich setzen. Weiter weg setzen
10. Welle mit einer ramp funktion anlaufen lassen
11. Diagonalen Messpunkt überprüfen und auf das gleiche delta setzen, wie den X Punkt. Wenn sich hier unterschiede auftun in der Messung bei längerem Abstand, stimmt die Messung nicht!
