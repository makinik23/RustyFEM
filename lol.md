# Kontekst projektu

Chcę zbudować edukacyjny, ale technicznie poprawny silnik metody elementów skończonych, dalej nazywany silnikiem MES lub FEM solverem, napisany w języku Rust.

Projekt ma służyć jednocześnie do:

1. praktycznej nauki języka Rust,
2. zrozumienia metody elementów skończonych od strony matematycznej,
3. poznania implementacji algorytmów numerycznych,
4. nauki projektowania oprogramowania inżynierskiego,
5. stworzenia wartościowego projektu do portfolio,
6. przygotowania fundamentu pod późniejsze zastosowania aerospace.

Docelowym zastosowaniem projektu mogą być między innymi:

* analiza prętów i kratownic,
* analiza cienkich paneli,
* analiza konstrukcji CubeSata,
* analiza paneli słonecznych,
* analiza wysięgników magnetometru,
* analiza płytek PCB,
* wyznaczanie przemieszczeń i naprężeń,
* analiza modalna i wyznaczanie częstości własnych.

Projekt nie ma być od początku odpowiednikiem komercyjnego solvera. Ma rozwijać się stopniowo, od bardzo prostych przypadków do bardziej zaawansowanych.

Najważniejsze są:

1. poprawność matematyczna,
2. poprawność numeryczna,
3. testowalność,
4. czytelność kodu,
5. zrozumiała architektura,
6. stopniowe rozszerzanie funkcjonalności,
7. dopiero na końcu wydajność i równoległość.

---

# Twoja rola

Działaj jednocześnie jako:

* senior Rust developer,
* inżynier metod numerycznych,
* architekt oprogramowania MES,
* mentor techniczny,
* recenzent kodu,
* konsultant mechaniki konstrukcji,
* autor testów i procedur walidacyjnych.

Nie ograniczaj się do generowania kodu.

Przed implementacją większego fragmentu:

1. opisz problem,
2. przedstaw potrzebną teorię,
3. zapisz równania,
4. określ wymiary macierzy i wektorów,
5. wyjaśnij kolejność stopni swobody,
6. zaproponuj strukturę danych,
7. zaproponuj API,
8. wskaż potencjalne błędy,
9. zaproponuj testy,
10. dopiero potem implementuj kod.

Jeżeli zauważysz błąd w moich założeniach matematycznych, fizycznych lub architektonicznych, nie realizuj go bezkrytycznie. Wskaż problem, wyjaśnij go i zaproponuj poprawne rozwiązanie.

---

# Główna filozofia projektu

Projekt ma być rozwijany w małych, działających przyrostach.

Każdy etap powinien kończyć się:

* działającym kodem,
* przechodzącymi testami,
* prostym przykładem użycia,
* dokumentacją założeń,
* walidacją wyniku,
* jasno opisanymi ograniczeniami.

Nie rozpoczynaj następnego etapu, dopóki obecny etap nie spełnia zdefiniowanych kryteriów ukończenia.

Nie próbuj przewidywać całej przyszłej architektury. Abstrakcje powinny wynikać z istniejących implementacji.

Nie stosuj nadmiernego overengineeringu.

---

# Docelowy zakres funkcjonalny

Docelowo solver powinien obsługiwać:

## Analizy

* liniową analizę statyczną,
* małe przemieszczenia,
* małe odkształcenia,
* analizę modalną.

## Materiały

* liniowy materiał sprężysty,
* moduł Younga,
* współczynnik Poissona,
* gęstość.

## Elementy

* pręt 1D,
* kratownicę 2D,
* element trójkątny T3/CST,
* element czworokątny Q4.

## Stany 2D

* płaski stan naprężenia,
* płaski stan odkształcenia.

## Obciążenia

* siły węzłowe,
* zadane przemieszczenia,
* obciążenia rozłożone na krawędziach,
* grawitację,
* opcjonalnie obciążenia termiczne w późniejszym etapie.

## Wyniki

* przemieszczenia,
* reakcje podporowe,
* odkształcenia,
* naprężenia,
* siły wewnętrzne,
* naprężenia von Misesa,
* częstości własne,
* postacie drgań.

## Algebra liniowa

* macierze gęste w pierwszych etapach,
* macierze rzadkie COO,
* macierze rzadkie CSR,
* bezpośredni solver układów liniowych,
* Conjugate Gradient,
* preconditioner Jacobiego.

## Wejście i wyjście

* prosty format modelu w JSON,
* import siatki z Gmsh,
* eksport wyników do VTK lub VTU,
* wizualizację wyników w ParaView,
* interfejs CLI.

---

# Funkcjonalności poza początkowym zakresem

Na początku nie implementuj:

* kontaktu,
* plastyczności,
* dużych odkształceń,
* nieliniowości materiałowej,
* nieliniowości geometrycznej,
* propagacji pęknięć,
* remeshingu,
* adaptacyjnej siatki,
* elementów bryłowych 3D,
* GPU,
* własnego GUI,
* własnego generatora siatki,
* zaawansowanych preconditionerów,
* metod wielosiatkowych,
* rozproszonego przetwarzania.

Każde z tych zagadnień może zostać dodane dopiero po ukończeniu stabilnej wersji podstawowej.

---

# Technologia

Projekt ma używać stabilnej wersji Rust.

Podstawowe narzędzia:

* Cargo,
* rustfmt,
* clippy,
* cargo test,
* cargo-nextest,
* cargo-llvm-cov,
* criterion do benchmarków.

Preferowane biblioteki:

* `nalgebra` do algebry liniowej,
* `thiserror` do obsługi błędów,
* `serde` do serializacji,
* `serde_json` do modeli JSON,
* `clap` do CLI,
* `rayon` dopiero na etapie równoległości,
* `approx` lub własna biblioteka tolerancji do testów numerycznych.

Nie implementuj kompletnej biblioteki algebry liniowej od zera jako podstawy finalnego solvera.

Możesz jednak edukacyjnie zaimplementować:

* prosty wektor,
* prostą macierz,
* eliminację Gaussa,
* rozkład LU,
* rozkład Cholesky'ego.

Po zakończeniu ćwiczeń właściwy solver powinien korzystać ze sprawdzonej biblioteki.

---

# Podstawowe zasady kodowania

Kod powinien być:

* idiomatyczny dla Rusta,
* bezpieczny,
* modularny,
* czytelny,
* testowalny,
* dobrze udokumentowany,
* pozbawiony nieuzasadnionych kopii,
* pozbawiony globalnego mutowalnego stanu.

Na początku używaj wyłącznie `f64`.

Nie generalizuj całego solvera na dowolny typ skalarny.

Nie używaj `unsafe`, chyba że na późniejszym etapie będzie to bardzo dobrze uzasadnione i poprzedzone benchmarkiem.

Nie używaj `unwrap()` ani `expect()` w kodzie biblioteki.

Publiczne funkcje powinny zwracać:

```rust
Result<T, FemError>
```

Dopuszczalne jest użycie `unwrap()` w:

* testach,
* małych przykładach,
* kodzie inicjalizującym CLI, jeżeli błąd jest wcześniej poprawnie opisany użytkownikowi.

Preferuj:

* małe funkcje,
* jawne typy domenowe,
* własne typy błędów,
* zwykłe pętle w kodzie numerycznym, gdy są bardziej czytelne niż złożone iteratory,
* iteratory tam, gdzie upraszczają kod,
* dokumentację `rustdoc`,
* testy jednostkowe blisko modułów,
* testy integracyjne w katalogu `tests`.

---

# Konwencje fizyczne i numeryczne

Wszystkie dane wejściowe i obliczenia powinny korzystać z układu SI.

Przyjmij:

* długość w metrach,
* pole w metrach kwadratowych,
* siłę w niutonach,
* moment w niutonometrach,
* moduł Younga w paskalach,
* naprężenie w paskalach,
* gęstość w kilogramach na metr sześcienny,
* masę w kilogramach,
* czas w sekundach,
* częstotliwość w hercach.

Indeksowanie wewnętrzne powinno zaczynać się od zera.

Identyfikatory w plikach wejściowych mogą mieć dowolne wartości i powinny być mapowane na indeksy wewnętrzne.

Dla problemów 2D przyjmij konwencję odkształceń:

[
\boldsymbol{\varepsilon}
========================

\begin{bmatrix}
\varepsilon_x \
\varepsilon_y \
\gamma_{xy}
\end{bmatrix},
]

gdzie (\gamma_{xy}) jest inżynierskim odkształceniem postaciowym.

Wektor naprężeń:

[
\boldsymbol{\sigma}
===================

\begin{bmatrix}
\sigma_x \
\sigma_y \
\tau_{xy}
\end{bmatrix}.
]

Konwencja ta musi być stosowana konsekwentnie w całym projekcie.

---

# Proponowana początkowa struktura repozytorium

Na początku projekt powinien być jednym crate'em zawierającym bibliotekę i program wykonywalny.

```text
rust-fem/
├── Cargo.toml
├── README.md
├── LICENSE
├── src/
│   ├── lib.rs
│   ├── main.rs
│   ├── error.rs
│   ├── math/
│   │   ├── mod.rs
│   │   └── tolerance.rs
│   ├── model/
│   │   ├── mod.rs
│   │   ├── node.rs
│   │   ├── material.rs
│   │   ├── section.rs
│   │   ├── load.rs
│   │   ├── constraint.rs
│   │   └── dof.rs
│   ├── elements/
│   │   ├── mod.rs
│   │   ├── bar_1d.rs
│   │   ├── truss_2d.rs
│   │   ├── triangle_t3.rs
│   │   └── quad_q4.rs
│   ├── assembly/
│   │   ├── mod.rs
│   │   ├── dense.rs
│   │   └── sparse.rs
│   ├── boundary/
│   │   ├── mod.rs
│   │   └── elimination.rs
│   ├── solver/
│   │   ├── mod.rs
│   │   ├── direct.rs
│   │   ├── conjugate_gradient.rs
│   │   └── preconditioner.rs
│   ├── analysis/
│   │   ├── mod.rs
│   │   ├── static_analysis.rs
│   │   └── modal_analysis.rs
│   ├── postprocessing/
│   │   ├── mod.rs
│   │   ├── stress.rs
│   │   ├── strain.rs
│   │   └── reactions.rs
│   └── io/
│       ├── mod.rs
│       ├── json.rs
│       ├── gmsh.rs
│       └── vtk.rs
├── tests/
│   ├── bar_1d.rs
│   ├── truss_2d.rs
│   ├── patch_test_t3.rs
│   ├── patch_test_q4.rs
│   └── modal_analysis.rs
├── examples/
├── models/
├── benches/
└── docs/
```

Nie twórz od razu workspace'u z wieloma crate'ami.

Podział na osobne crate'y może nastąpić dopiero wtedy, gdy:

* moduły będą stabilne,
* zależności między nimi będą dobrze zrozumiane,
* publiczne API będzie dojrzałe,
* jeden crate rzeczywiście zacznie utrudniać rozwój.

---

# Obsługa błędów

Zaprojektuj własny typ błędu.

Przykładowo:

```rust
#[derive(Debug, thiserror::Error)]
pub enum FemError {
    #[error("node with id {0} does not exist")]
    InvalidNodeId(usize),

    #[error("material with id {0} does not exist")]
    InvalidMaterialId(usize),

    #[error("section with id {0} does not exist")]
    InvalidSectionId(usize),

    #[error("element with id {0} does not exist")]
    InvalidElementId(usize),

    #[error("invalid element geometry: {0}")]
    InvalidElementGeometry(String),

    #[error("matrix dimensions are inconsistent")]
    InvalidMatrixDimensions,

    #[error("global stiffness matrix is singular")]
    SingularMatrix,

    #[error("model contains unconstrained rigid body modes")]
    UnconstrainedRigidBodyModes,

    #[error("linear solver did not converge")]
    SolverDidNotConverge,

    #[error("invalid material properties: {0}")]
    InvalidMaterialProperties(String),

    #[error("input data error: {0}")]
    InputData(String),
}
```

Błędy powinny zawierać wystarczająco dużo informacji, aby użytkownik mógł zlokalizować problem.

Przykładowo błąd odwróconego elementu Q4 powinien podawać:

* identyfikator elementu,
* punkt całkowania,
* wartość wyznacznika Jacobianu.

---

# Etap 0 — inicjalizacja projektu

Zadania:

1. utworzenie projektu Cargo,
2. utworzenie biblioteki i programu CLI,
3. konfiguracja `rustfmt`,
4. konfiguracja `clippy`,
5. konfiguracja testów,
6. utworzenie podstawowych modułów,
7. przygotowanie README,
8. konfiguracja CI,
9. określenie zasad jednostek i indeksowania.

Minimalne zależności mogą obejmować:

```toml
[dependencies]
nalgebra = "..."
thiserror = "..."
serde = { version = "...", features = ["derive"] }
serde_json = "..."
clap = { version = "...", features = ["derive"] }

[dev-dependencies]
approx = "..."
criterion = "..."
```

Nie wpisuj numerów wersji z pamięci. Przed dodaniem zależności sprawdź aktualne, stabilne wersje zgodne z używaną wersją Rust.

Kryteria ukończenia:

```bash
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test
```

Wszystkie polecenia powinny przechodzić bez błędów.

---

# Etap 1 — edukacyjna algebra liniowa

Ten etap ma służyć nauce Rusta i metod numerycznych.

Zaimplementuj proste typy:

```rust
pub struct Vector {
    data: Vec<f64>,
}

pub struct Matrix {
    rows: usize,
    cols: usize,
    data: Vec<f64>,
}
```

Dodaj:

* konstruktory,
* kontrolę wymiarów,
* bezpieczne indeksowanie,
* dodawanie,
* odejmowanie,
* mnożenie przez skalar,
* mnożenie macierz–wektor,
* mnożenie macierz–macierz,
* transpozycję,
* iloczyn skalarny,
* normę euklidesową.

Następnie edukacyjnie zaimplementuj:

* eliminację Gaussa z pivotingiem,
* rozkład LU,
* forward substitution,
* backward substitution,
* rozkład Cholesky'ego.

Dla każdego algorytmu:

* opisz założenia,
* podaj złożoność czasową,
* opisz stabilność numeryczną,
* wykrywaj macierze osobliwe,
* przygotuj testy.

Po ukończeniu etapu kod produkcyjny solvera powinien przejść na `nalgebra`.

Nie rozwijaj własnej algebry liniowej w kompletny zamiennik `nalgebra`.

---

# Etap 2 — pręt 1D

Pierwszy pełny solver ma rozwiązywać liniowy problem osiowo obciążonego pręta.

Równanie silne:

[
\frac{d}{dx}
\left(
EA\frac{du}{dx}
\right)
+
f(x)
====

0.

]

Dla dwuwęzłowego elementu o stałych parametrach:

[
\mathbf{k}^{(e)}
================

\frac{EA}{L}
\begin{bmatrix}
1 & -1 \
-1 & 1
\end{bmatrix}.
]

Struktury domenowe powinny rozdzielać:

* węzeł,
* materiał,
* przekrój,
* element,
* obciążenie,
* więzy.

Przykład:

```rust
pub struct Node1D {
    id: usize,
    x: f64,
}

pub struct LinearElasticMaterial {
    id: usize,
    young_modulus: f64,
    poisson_ratio: f64,
    density: f64,
}

pub struct BarSection {
    id: usize,
    area: f64,
}

pub struct BarElement1D {
    id: usize,
    node_ids: [usize; 2],
    material_id: usize,
    section_id: usize,
}
```

Zaimplementuj:

1. walidację parametrów materiału,
2. walidację pola przekroju,
3. obliczenie długości elementu,
4. wykrywanie zerowej długości,
5. lokalną macierz sztywności,
6. mapowanie lokalnych DOF na globalne,
7. assembly macierzy globalnej,
8. assembly globalnego wektora sił,
9. zadane przemieszczenia,
10. rozwiązanie zredukowanego układu,
11. odtworzenie pełnego wektora przemieszczeń,
12. reakcje podporowe,
13. odkształcenia,
14. naprężenia,
15. siły osiowe.

Dla elementu:

[
\varepsilon
===========

\frac{u_2-u_1}{L},
]

[
\sigma
======

E\varepsilon,
]

[
N
=

A\sigma.
]

Podstawowy przypadek testowy:

* pręt długości (L),
* pole przekroju (A),
* moduł Younga (E),
* lewy koniec utwierdzony,
* siła (F) na prawym końcu.

Oczekiwane przemieszczenie:

[
u(L)
====

\frac{FL}{EA}.
]

Dodatkowe testy:

* dwa elementy o tych samych właściwościach,
* dwa elementy o różnych polach,
* dwa elementy o różnych materiałach,
* niezerowe przemieszczenie zadane,
* obciążenie pośredniego węzła,
* brak wystarczających więzów,
* element o zerowej długości,
* niepoprawny identyfikator materiału,
* równowaga reakcji i sił.

Kryterium równowagi:

[
\sum F_i+\sum R_i \approx 0.
]

---

# Etap 3 — uporządkowanie architektury domenowej

Dopiero po działającym solverze pręta 1D zaproponuj refactor.

Możliwe typy:

* `Node`,
* `Material`,
* `Section`,
* `Element`,
* `Load`,
* `Constraint`,
* `Model`,
* `Dof`,
* `Solution`,
* `StaticAnalysisResult`.

Przykładowy trait:

```rust
pub trait FiniteElement {
    fn node_ids(&self) -> &[usize];

    fn dof_indices(
        &self,
        model: &Model,
    ) -> Result<Vec<usize>, FemError>;

    fn stiffness_matrix(
        &self,
        model: &Model,
    ) -> Result<ElementMatrix, FemError>;
}
```

Nie zakładaj, że wszystkie przyszłe elementy będą miały identyczne potrzeby.

Trait powinien być możliwie mały.

Rozważ także użycie `enum` dla typów elementów, jeżeli będzie to prostsze i łatwiejsze do serializacji niż dynamiczny dispatch.

Przed wyborem między:

* `enum Element`,
* `Box<dyn FiniteElement>`,
* typami generycznymi,

porównaj zalety i wady.

Dla początkowego solvera rekomenduj rozwiązanie najprostsze.

---

# Etap 4 — kratownica 2D

Każdy węzeł ma dwa stopnie swobody:

[
u_x,\quad u_y.
]

Kolejność lokalnych DOF:

[
\mathbf u_e
===========

\begin{bmatrix}
u_{1x} &
u_{1y} &
u_{2x} &
u_{2y}
\end{bmatrix}^T.
]

Dla elementu:

[
c
=

\frac{x_2-x_1}{L},
\qquad
s
=

\frac{y_2-y_1}{L}.
]

Globalna macierz sztywności elementu:

[
\mathbf{k}^{(e)}
================

\frac{EA}{L}
\begin{bmatrix}
c^2 & cs & -c^2 & -cs \
cs & s^2 & -cs & -s^2 \
-c^2 & -cs & c^2 & cs \
-cs & -s^2 & cs & s^2
\end{bmatrix}.
]

Zaimplementuj:

* węzły 2D,
* mapowanie dwóch DOF na węzeł,
* element kratownicowy,
* siły osiowe,
* naprężenia,
* reakcje,
* kontrolę mechanizmów.

Testy:

* element poziomy,
* element pionowy,
* element ukośny,
* kratownica trójkątna,
* model statycznie wyznaczalny,
* model z mechanizmem,
* test obrotu całej geometrii,
* test translacji całej geometrii.

Po obróceniu geometrii, obciążeń i więzów rozwiązanie powinno obrócić się zgodnie z tą samą transformacją.

---

# Etap 5 — element T3/CST

Zaimplementuj liniowy element trójkątny dla:

* płaskiego stanu naprężenia,
* płaskiego stanu odkształcenia.

Macierz sztywności:

[
\mathbf K^{(e)}
===============

tA\mathbf B^T\mathbf D\mathbf B.
]

Dla płaskiego stanu naprężenia:

[
\mathbf D
=========

\frac{E}{1-\nu^2}
\begin{bmatrix}
1 & \nu & 0 \
\nu & 1 & 0 \
0 & 0 & \frac{1-\nu}{2}
\end{bmatrix}.
]

Dla płaskiego stanu odkształcenia:

[
\mathbf D
=========

\frac{E}{(1+\nu)(1-2\nu)}
\begin{bmatrix}
1-\nu & \nu & 0 \
\nu & 1-\nu & 0 \
0 & 0 & \frac{1-2\nu}{2}
\end{bmatrix}.
]

Zaimplementuj:

* pole elementu,
* orientację węzłów,
* funkcje kształtu,
* pochodne funkcji kształtu,
* macierz (\mathbf B),
* macierz konstytutywną (\mathbf D),
* macierz sztywności,
* odkształcenia,
* naprężenia,
* von Mises.

Odkształcenia:

[
\boldsymbol{\varepsilon}
========================

\mathbf B\mathbf u_e.
]

Naprężenia:

[
\boldsymbol{\sigma}
===================

\mathbf D\boldsymbol{\varepsilon}.
]

Dla płaskiego stanu naprężenia:

[
\sigma_{vm}
===========

\sqrt{
\sigma_x^2
----------

\sigma_x\sigma_y
+
\sigma_y^2
+
3\tau_{xy}^2
}.
]

Obowiązkowe testy:

* poprawność pola,
* suma funkcji kształtu,
* kompletność liniowa,
* symetria macierzy sztywności,
* ruch ciała sztywnego,
* patch test,
* jednorodne rozciąganie,
* odwrócona kolejność węzłów,
* element zdegenerowany.

Patch test jest warunkiem przejścia do kolejnego etapu.

---

# Etap 6 — element Q4

Zaimplementuj czterowęzłowy element izoparametryczny Q4.

Funkcje kształtu:

[
N_1
===

\frac14(1-\xi)(1-\eta),
]

[
N_2
===

\frac14(1+\xi)(1-\eta),
]

[
N_3
===

\frac14(1+\xi)(1+\eta),
]

[
N_4
===

\frac14(1-\xi)(1+\eta).
]

Zaimplementuj:

* pochodne po (\xi) i (\eta),
* Jacobian,
* wyznacznik Jacobianu,
* odwrotność Jacobianu,
* transformację pochodnych do układu globalnego,
* macierz (\mathbf B),
* kwadraturę Gaussa (2\times2),
* macierz sztywności,
* naprężenia w punktach Gaussa.

W każdym punkcie całkowania wymagaj:

[
\det(\mathbf J)>0.
]

Jeżeli warunek nie jest spełniony, zwróć opisowy błąd.

Testy:

* suma funkcji kształtu,
* kompletność liniowa,
* regularny prostokąt,
* zniekształcony czworokąt,
* patch test,
* odwrócona numeracja,
* element zdegenerowany,
* porównanie wyników Q4 i T3.

---

# Etap 7 — warunki brzegowe

Domyślną metodą nakładania zadanych przemieszczeń ma być redukcja układu.

Podział:

[
\begin{bmatrix}
K_{ff} & K_{fc} \
K_{cf} & K_{cc}
\end{bmatrix}
\begin{bmatrix}
u_f \
u_c
\end{bmatrix}
=============

\begin{bmatrix}
F_f \
F_c
\end{bmatrix}.
]

Układ dla swobodnych DOF:

[
K_{ff}u_f
=========

F_f-K_{fc}u_c.
]

Obsłuż:

* zerowe zadane przemieszczenia,
* niezerowe zadane przemieszczenia,
* wiele więzów,
* konfliktujące więzy,
* więzy na nieistniejącym DOF.

Reakcje wyznaczaj z pełnego, oryginalnego układu:

[
\mathbf R
=========

\mathbf K\mathbf u-\mathbf F.
]

Nie obliczaj reakcji z macierzy zmodyfikowanej przez więzy.

---

# Etap 8 — obciążenia

Obsłuż:

* siły węzłowe,
* siły rozłożone na krawędzi,
* grawitację,
* opcjonalnie obciążenia termiczne.

Dla obciążenia krawędziowego:

[
\mathbf f_e
===========

\int_{\Gamma_e}
\mathbf N^T\mathbf t,d\Gamma.
]

Dla każdego rodzaju obciążenia przygotuj:

* model danych,
* walidację,
* równoważny wektor sił,
* test analityczny,
* test równowagi.

---

# Etap 9 — macierze rzadkie

Zachowaj wersję gęstą jako implementację referencyjną.

Dodaj:

* format COO,
* sumowanie duplikatów,
* sortowanie wpisów,
* konwersję COO do CSR,
* mnożenie CSR przez wektor,
* sparse assembly.

Przykładowe struktury:

```rust
pub struct CooMatrix {
    nrows: usize,
    ncols: usize,
    rows: Vec<usize>,
    cols: Vec<usize>,
    values: Vec<f64>,
}

pub struct CsrMatrix {
    nrows: usize,
    ncols: usize,
    row_offsets: Vec<usize>,
    column_indices: Vec<usize>,
    values: Vec<f64>,
}
```

Testy:

* porównanie z macierzą gęstą,
* sumowanie duplikatów,
* puste wiersze,
* mnożenie przez wektor,
* assembly tego samego modelu przez dense i sparse,
* zgodność rozwiązania.

---

# Etap 10 — Conjugate Gradient

Zaprojektuj interfejs operatora liniowego:

```rust
pub trait LinearOperator {
    fn dimension(&self) -> usize;

    fn apply(
        &self,
        x: &[f64],
        y: &mut [f64],
    ) -> Result<(), FemError>;
}
```

Zaimplementuj:

* Conjugate Gradient,
* preconditioner Jacobiego,
* kryterium zbieżności,
* limit iteracji,
* raport residualu,
* wykrywanie stagnacji.

Residual:

[
\mathbf r_k
===========

\mathbf b-\mathbf A\mathbf x_k.
]

Kryterium względne:

[
\frac{
|\mathbf r_k|_2
}{
\max(|\mathbf b|_2,\varepsilon)
}
<
\text{tolerance}.
]

Wynik solvera:

```rust
pub struct SolverResult {
    pub solution: Vec<f64>,
    pub converged: bool,
    pub iterations: usize,
    pub residual_norm: f64,
    pub relative_residual_norm: f64,
}
```

Obsłuż:

* zerowy wektor prawej strony,
* zbieżność początkowego przybliżenia,
* brak dodatniej określoności,
* dzielenie przez wartość bliską zeru,
* brak zbieżności,
* osiągnięcie limitu iteracji.

Porównaj:

* solver bezpośredni,
* CG,
* preconditioned CG.

---

# Etap 11 — format JSON i CLI

Model wejściowy powinien zawierać:

* typ analizy,
* konfigurację solvera,
* węzły,
* materiały,
* przekroje,
* elementy,
* więzy,
* obciążenia.

Przykład:

```json
{
  "analysis": {
    "type": "linear_static",
    "solver": "direct"
  },
  "nodes": [
    {
      "id": 1,
      "coordinates": [0.0]
    },
    {
      "id": 2,
      "coordinates": [1.0]
    }
  ],
  "materials": [
    {
      "id": 1,
      "type": "linear_elastic",
      "young_modulus": 70000000000.0,
      "poisson_ratio": 0.33,
      "density": 2700.0
    }
  ],
  "sections": [
    {
      "id": 1,
      "type": "bar",
      "area": 0.0001
    }
  ],
  "elements": [
    {
      "id": 1,
      "type": "bar_1d",
      "nodes": [1, 2],
      "material": 1,
      "section": 1
    }
  ],
  "constraints": [
    {
      "node": 1,
      "dof": "ux",
      "value": 0.0
    }
  ],
  "loads": [
    {
      "node": 2,
      "dof": "ux",
      "value": 1000.0
    }
  ]
}
```

CLI powinno umożliwiać:

```bash
cargo run --release -- solve models/bar.json
```

Raport powinien zawierać:

* liczbę węzłów,
* liczbę elementów,
* liczbę DOF,
* liczbę więzów,
* typ solvera,
* czas assembly,
* czas rozwiązania,
* residual,
* maksymalne przemieszczenie,
* lokalizację plików wynikowych.

---

# Etap 12 — Gmsh i VTK

Nie implementuj własnego meshera.

Importuj siatki z Gmsh.

Na początku obsłuż:

* węzły,
* linie,
* trójkąty,
* czworokąty,
* physical groups.

Eksportuj do VTK lub VTU:

## Dane węzłowe

* przemieszczenia,
* wartość przemieszczenia,
* reakcje.

## Dane elementowe

* odkształcenia,
* naprężenia,
* von Mises,
* identyfikator elementu,
* identyfikator materiału.

Eksport powinien zachować:

* geometrię niezdeformowaną,
* wektor przemieszczenia,
* możliwość wizualizacji deformacji w ParaView.

---

# Etap 13 — analiza modalna

Rozwiązuj uogólniony problem własny:

[
\mathbf K\boldsymbol{\phi}
==========================

\omega^2\mathbf M\boldsymbol{\phi}.
]

Dodaj:

* macierz masy,
* masę skupioną,
* masę konsystentną,
* assembly macierzy masy,
* redukcję warunków brzegowych,
* rozwiązanie problemu własnego,
* sortowanie częstotliwości,
* normalizację postaci drgań,
* eksport modów.

Częstotliwość:

[
f_i
===

\frac{\omega_i}{2\pi}.
]

Dla każdego modu oblicz residual:

[
\mathbf r_i
===========

## \mathbf K\boldsymbol{\phi}_i

\omega_i^2\mathbf M\boldsymbol{\phi}_i.
]

Sprawdzaj ortogonalność względem macierzy masy:

[
\boldsymbol{\phi}_i^T
\mathbf M
\boldsymbol{\phi}_j
\approx 0
\quad
\text{dla } i\neq j.
]

Waliduj wyniki na:

* pręcie osiowym,
* prostej kratownicy,
* prostym panelu,
* modelu referencyjnym z MATLAB-a lub innego solvera.

---

# Strategia testowania

Testy są częścią implementacji, a nie zadaniem wykonywanym na końcu.

## Testy jednostkowe

Testuj:

* konstruktory,
* walidację,
* funkcje kształtu,
* pochodne funkcji kształtu,
* Jacobian,
* macierze materiałowe,
* macierze elementowe,
* mapowanie DOF,
* assembly,
* warunki brzegowe,
* solvery,
* macierze sparse.

## Testy matematycznych własności

Sprawdzaj:

[
\mathbf K
=========

\mathbf K^T.
]

Sprawdzaj energię:

[
\mathbf u^T\mathbf K\mathbf u
\geq 0
]

z uwzględnieniem tolerancji.

Dla funkcji kształtu:

[
\sum_i N_i
==========

1.

]

Dla ich pochodnych:

[
\sum_i
\frac{\partial N_i}{\partial x}
===============================

0,
]

[
\sum_i
\frac{\partial N_i}{\partial y}
===============================

0.

]

Sprawdzaj także:

* niezmienniczość względem translacji geometrii,
* poprawną transformację względem obrotu,
* zerową energię dla ruchów ciała sztywnego przed nałożeniem więzów,
* równowagę sił i reakcji.

## Testy integracyjne

Modele referencyjne:

* pojedynczy pręt,
* szereg prętów,
* kratownica trójkątna,
* patch test T3,
* patch test Q4,
* rozciągana płyta,
* belka wspornikowa,
* panel do analizy modalnej.

## Testy regresyjne

Zapisuj i porównuj:

* wybrane przemieszczenia,
* reakcje,
* naprężenia,
* residual,
* liczbę iteracji,
* częstotliwości własne.

---

# Porównywanie liczb zmiennoprzecinkowych

Nie używaj:

```rust
assert_eq!(a, b);
```

dla wyników numerycznych.

Stosuj:

[
|a-b|
\leq
\varepsilon_{abs}
+
\varepsilon_{rel}
\max(|a|,|b|).
]

Przygotuj wspólną funkcję:

```rust
pub fn approx_eq(
    actual: f64,
    expected: f64,
    abs_tol: f64,
    rel_tol: f64,
) -> bool;
```

Tolerancje powinny być dobierane świadomie do konkretnego testu.

Nie zwiększaj tolerancji tylko po to, aby ukryć błąd implementacji.

---

# Dokumentacja

Każdy moduł elementu powinien dokumentować:

* fizyczne znaczenie elementu,
* założenia,
* kolejność węzłów,
* kolejność lokalnych DOF,
* układ współrzędnych,
* równania,
* jednostki,
* ograniczenia,
* przypadki błędne,
* test referencyjny.

Przykładowo dla Q4 dokumentacja powinna jednoznacznie wskazywać oczekiwaną kolejność węzłów.

Każda publiczna funkcja powinna mieć dokumentację `rustdoc`.

---

# Profilowanie i optymalizacja

Nie optymalizuj przed uzyskaniem poprawnych wyników.

Gdy solver będzie poprawny, mierz osobno:

* czas wczytania modelu,
* czas budowy mapowania DOF,
* czas assembly,
* czas nakładania więzów,
* czas rozwiązania,
* czas postprocessingu,
* czas eksportu.

Dopiero na podstawie pomiarów rozważ:

* prealokację,
* redukcję liczby alokacji,
* użycie statycznych macierzy dla małych elementów,
* zmianę formatu sparse,
* równoległy postprocessing,
* równoległe assembly,
* `rayon`,
* bardziej zaawansowane preconditionery.

Nie zmieniaj kodu na mniej czytelny bez wykazanego zysku wydajności.

---

# Praca z istniejącym repozytorium

Przed każdą większą zmianą:

1. przeanalizuj strukturę repozytorium,
2. przeczytaj `Cargo.toml`,
3. przeczytaj README,
4. sprawdź publiczne API,
5. przeczytaj istniejące testy,
6. uruchom testy,
7. sprawdź historię i styl istniejącego kodu,
8. unikaj przepisywania działających fragmentów bez powodu.

Jeżeli konieczny jest refactor:

* wyjaśnij jego cel,
* oddziel go od nowej funkcjonalności,
* nie zmieniaj zachowania bez testów,
* zachowaj kompatybilność API, jeżeli jest to rozsądne.

---

# Workflow dla każdego zadania

Dla każdego zadania stosuj poniższy schemat.

## 1. Analiza

Najpierw podaj:

* cel zadania,
* aktualny stan projektu,
* zależności,
* założenia,
* ryzyka.

## 2. Matematyka

Podaj:

* równania,
* wymiary macierzy,
* kolejność DOF,
* oczekiwane własności wyników.

## 3. Projekt rozwiązania

Zaproponuj:

* typy,
* funkcje,
* interfejsy,
* moduły,
* typy błędów.

## 4. Plan testów

Przed implementacją opisz:

* test jednostkowy,
* test analityczny,
* test błędnych danych,
* test integracyjny,
* kryterium akceptacji.

## 5. Implementacja

Wprowadzaj możliwie mały przyrost.

Nie zmieniaj plików niezwiązanych z zadaniem.

Nie generuj dużej liczby abstrakcji przed ich wykorzystaniem.

## 6. Weryfikacja

Uruchom:

```bash
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test
```

W razie potrzeby:

```bash
cargo nextest run
cargo llvm-cov
cargo bench
```

## 7. Podsumowanie

Na końcu każdej zmiany podaj:

* listę zmienionych plików,
* opis działania,
* przyjęte decyzje,
* uruchomione testy,
* wyniki testów,
* znane ograniczenia,
* następny najmniejszy sensowny krok.

---

# Zasady współpracy

Nie implementuj całego roadmapu jednym ogromnym commitem.

Nie generuj kilkudziesięciu plików bez wcześniejszego uzgodnienia architektury.

Nie twórz placeholderów typu `todo!()` w dużej liczbie modułów, które nie są jeszcze potrzebne.

Nie implementuj „na zapas” funkcjonalności, które będą potrzebne dopiero za wiele etapów.

Jeżeli istnieje kilka możliwych rozwiązań:

1. przedstaw maksymalnie trzy,
2. porównaj zalety i wady,
3. wybierz rozwiązanie rekomendowane,
4. uzasadnij wybór w kontekście aktualnego etapu.

Preferuj rozwiązanie najprostsze, które:

* jest poprawne,
* jest testowalne,
* pozwala przejść do następnego etapu,
* nie zamyka oczywistej drogi rozwoju.

---

# Definition of Done dla nowego elementu MES

Nowy typ elementu można uznać za ukończony dopiero wtedy, gdy posiada:

1. opis matematyczny,
2. określoną kolejność węzłów,
3. określoną kolejność DOF,
4. walidację geometrii,
5. lokalną macierz sztywności,
6. mapowanie DOF,
7. postprocessing,
8. test jednostkowy,
9. test rozwiązania analitycznego lub referencyjnego,
10. patch test, jeżeli ma zastosowanie,
11. test błędnej geometrii,
12. dokumentację,
13. przykład użycia.

---

# Pierwszy milestone

Pierwszym rzeczywistym milestone'em ma być:

## Liniowy solver statyczny dla prętów 1D

Zakres milestone'u:

* dowolna liczba węzłów,
* dowolna liczba elementów,
* wiele materiałów,
* wiele przekrojów,
* siły węzłowe,
* zerowe zadane przemieszczenia,
* niezerowe zadane przemieszczenia,
* globalny assembly,
* redukcja układu,
* rozwiązanie,
* rekonstrukcja pełnego wektora przemieszczeń,
* reakcje podporowe,
* odkształcenia,
* naprężenia,
* siły osiowe,
* testy analityczne,
* prosty model JSON,
* minimalny interfejs CLI.

Na tym etapie nie implementuj jeszcze:

* kratownic 2D,
* macierzy rzadkich,
* VTK,
* Gmsh,
* analizy modalnej.

---

# Pierwsze zadanie dla Codexa

Rozpocznij od analizy i projektu pierwszego milestone'u.

W pierwszej odpowiedzi:

1. przeanalizuj obecne repozytorium, jeżeli już istnieje,
2. zaproponuj minimalną strukturę projektu,
3. zaproponuj zależności w `Cargo.toml`,
4. zdefiniuj minimalne typy domenowe,
5. opisz przepływ danych od modelu wejściowego do wyniku,
6. przedstaw lokalną macierz sztywności pręta,
7. opisz mapowanie lokalnych DOF na globalne,
8. opisz assembly,
9. opisz nakładanie więzów metodą redukcji,
10. opisz obliczanie reakcji,
11. opisz postprocessing elementu,
12. zaproponuj testy,
13. podziel milestone na małe, kolejno wykonywane zadania,
14. wskaż ryzyka architektoniczne i numeryczne.

Nie implementuj od razu całego solvera.

Po przedstawieniu planu rozpocznij implementację wyłącznie pierwszego, najmniejszego kroku, którym powinno być przygotowanie podstawowych typów domenowych i ich walidacji wraz z testami.

Po każdym kroku zachowaj działające:

```bash
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test
```
