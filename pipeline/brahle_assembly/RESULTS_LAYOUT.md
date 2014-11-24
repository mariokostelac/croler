Rezulati pokretanja
===================

Nacin na koji sam testirao je:

* Generirao sam readove koristeci http://sourceforge.net/projects/readsim/
* Koristeci amos, odnosno minimus, stvorio sam reads.afg i overlaps.afg datoteke
* Pokrenuo sam `bin/main_layout` programa na tako stvorenim ulazima
* Vrijeme za pojedine faze program ispisuje sam, unitigging faza (glavni dio algoritma), prikazana je tablici u stupcu "Samo unitigging"
* "Kompletno rjesenje" je "stvarno" (real) vrijeme koje sam dobio mjereci duzinu izvrsavanja programa koristeci unix naredbu `/usr/bin/time`

Rezultati
-------

| Primjer | Readova | Overlapa | Minimus (faza_layouta) | Samo unitigging | Kompletno rjesenje | n50 | Najveci kontig
|---------|---------|----------|------------------------|-----------------|--------------------|-----|----------------
| 1       |     167 |     1323 |                    0 s |          0.00 s |             0.03 s |  21 |             21 
| 2       |   23856 |   249674 |                    2 s |          0.09 s |             3.53 s | 101 |            225


Komentari o primjerima
----------------------

* Primjer 1
 * Prvih 500 redaka iz readsim-1.5/example/ecoli/NC_000913.fna, sto je ukupno 34390 baza. 
 * I minimus i moje rjesenje dobije isti kontig
* Primjer 2
 * Cijela readsim-1.5/example/ecoli/NC_000913.fna datoteka, sto je ukupno 4638675 baza.
 * Minimus pronadje 47 kontiga, ovo rjesnje pronadje 43
  * Razlog je vjerojatno sto ovo rjesenje agresivnije brise tranzitivne bridove. 
 * Glavnina vremena potrosena je na IO

Opcenito komentari
------------------

Koristio sam callgrind (`valgrind --tool=callgrind bin/main_layout` i onda sam pregledavao datoteku koristeci `kcachegrind`) da vidim zasto kompletano rjesenje traje toliko puno. Najsporiji dio u rjesenju je IO - na to se potrosi 90% vremena, prema callgrindu. Kada spojimo rjesenje, problem s IO ce nestati, tako da bi to trebala biti najveca optimizacija. Sve 3 faze unitigginga su otprilike jednako brze.
