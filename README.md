OpenGL-DancingControlPoints
===========================

OpenGL homework with curves

Feladat
------------
Készítsen *„Síkon táncoló kontrollpontok”* programot.
A felhasználó az egér *balgomb* egyszeri lenyomással/elengedéssel veszi fel a kontrollpontokat (max 10-et),
amelyekhez 2cm sugarú kisebb fekete köröket rendelünk.

Ha a kontrollpontok száma legalább kettő,
azokra egy türkiszkék színű, kitöltött **konvex burkot**,

piros **Bézier görbét**,

nulla kezdő és végsebességű, a kontrollpont lehelyezésének idejét paraméterként használó zöld **Catmull-Rom spline-t**

és kék **Catmull-Clark görbét** illeszt.


A háttér világosszürke.
Legnagyobb prioritása a kontrollpontoknak van, majd a görbék jönnek, végül jön a konvex burok.

*Space* lenyomására a kontrollpontok egy-egy 5 cm-es kör tetejéről elindulva,
5 másodpercenként egy teljes fordulatot téve, elkezdenek keringeni,
mégpedig a páros indexűek az óramutató járásával megegyező, a páratlan indexűek pedig azzal ellentétes irányban.

A konvex burok és a görbék követik a kontrollpontokat.
Mindezt a felhasználó egy 58cm x 68 cm-es kameraablakon keresztül látja.

Ha a felhasználó az egér *jobb gombbal* rábök egy kontrollpontra, akkor a kameraablakot ehhez köti,
a kontrollpont elmozdulása automatikusan a kameraablakot is arrébb viszi.

Beadási határidő: 2014. 10. 19. 23:59
