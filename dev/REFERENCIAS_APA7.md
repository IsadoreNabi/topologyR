# Las referencias del paquete, en APA 7: lista canónica

Cada entrada se escribe una vez acá y se copia literal a los bloques `@references` de Roxygen; un bloque que no
coincida con esta lista es un defecto, y `dev/verificar_referencias.R` lo comprueba. La carpeta `dev/` está
excluida de la construcción del paquete por `.Rbuildignore`.

## Decisiones de forma

1. Título en caja de oración (la primera palabra, la primera después de dos puntos y los nombres propios), como
   pide APA 7 para artículos y libros; el sustantivo alemán de Alexandroff conserva su mayúscula.
2. Los localizadores van en la prosa de `@details` o de las secciones, nunca en la referencia.
3. DOI en forma de dirección `https://doi.org/...` cuando existe; sin DOI, la dirección del archivo oficial.
4. Los nombres van con sus diacríticos (Nuño, Nuñez, Roldán): el DESCRIPTION declara `Encoding: UTF-8` y los
   diacríticos sólo aparecen en comentarios de Roxygen y en los `.Rd` generados, no en el código.

## Verificación de los datos bibliográficos

- Consulta a la interfaz de Crossref del 2026-10-05 a las 09:20:14 -06:00 (`refs/crossref_dois.txt` del registro de
  R02): Nada y otros (2018), Lacasa y otros (2008), Luque y otros (2009), Lacasa y otros (2012), Kelly (1963), Tarjan
  (1972), Ganter y Wille (1999) y Shewchuk (1997); autores, revista, volumen, número y páginas coinciden con lo de abajo.
- Alexandroff (1937): ficha del archivo oficial Math-Net.Ru, https://www.mathnet.ru/eng/sm5579, consultada el
  2026-10-05: «Rec. Math. [Mat. Sbornik] N.S.», 1937, 2(44), número 3, 501–519. Un resultado de buscador daba
  501–518; manda el archivo oficial.
- Las obras depositadas en `bitopology-nada/BIBLIOGRAFÍA` (Nada y otros, 2018; Lacasa y otros, 2012) se cotejaron
  además con su propio PDF.

## Artículos de revista (§10.1)

    Alexandroff, P. (1937). Diskrete Räume. Recueil Mathématique (Matematicheskii Sbornik), Nouvelle Série,
      2(44)(3), 501-519. https://www.mathnet.ru/eng/sm5579

    Kelly, J. C. (1963). Bitopological spaces. Proceedings of the London Mathematical Society, s3-13(1), 71-89.
      https://doi.org/10.1112/plms/s3-13.1.71

    Lacasa, L., Luque, B., Ballesteros, F., Luque, J., & Nuño, J. C. (2008). From time series to complex networks:
      The visibility graph. Proceedings of the National Academy of Sciences, 105(13), 4972-4975.
      https://doi.org/10.1073/pnas.0709247105

    Lacasa, L., Nuñez, A., Roldán, É., Parrondo, J. M. R., & Luque, B. (2012). Time series irreversibility: A
      visibility graph approach. The European Physical Journal B, 85(6), Article 217.
      https://doi.org/10.1140/epjb/e2012-20809-8

    Luque, B., Lacasa, L., Ballesteros, F., & Luque, J. (2009). Horizontal visibility graphs: Exact results for
      random time series. Physical Review E, 80(4), Article 046103. https://doi.org/10.1103/PhysRevE.80.046103

    Nada, S., El Atik, A. E. F., & Atef, M. (2018). New types of topological structures via graphs. Mathematical
      Methods in the Applied Sciences, 41(15), 5801-5810. https://doi.org/10.1002/mma.4726

    Shewchuk, J. R. (1997). Adaptive precision floating-point arithmetic and fast robust geometric predicates.
      Discrete & Computational Geometry, 18(3), 305-363. https://doi.org/10.1007/PL00009321

    Tarjan, R. (1972). Depth-first search and linear graph algorithms. SIAM Journal on Computing, 1(2), 146-160.
      https://doi.org/10.1137/0201010

## Libros (§10.2)

    Ganter, B., & Wille, R. (1999). Formal concept analysis: Mathematical foundations. Springer.
      https://doi.org/10.1007/978-3-642-59830-2

## Referencias retiradas en 0.4.0, y por qué

- Lacasa y Toral (2010): estaba citada como fuente de una «predicción de irreversibilidad»; su texto completo no
  contiene «irrevers» (lectura de la ronda 7). El trabajo pertinente sobre irreversibilidad con grafos de
  visibilidad es Lacasa y otros (2012).
- Birkhoff (1940): estaba citada como origen de la dualidad de Galois que usa la prueba del teorema de los
  cardinales de base. La prueba queda escrita completa en la documentación y la atribución histórica no se
  verificó contra la obra, así que no se afirma.
- La atribución de la conexidad por pares a Kelly (1963): los dos textos leídos que la definen (Abdu y Kılıçman,
  2018; Baby Girija y Pilakkat, 2013) no la atribuyen a Kelly, de modo que la documentación la define sin
  atribuirla, y cita a Kelly sólo como origen de los espacios bitopológicos, que los dos textos le atribuyen.
