# kodesChemistryModel

Библиотека OpenFOAM для расчёта химических реакций на GPU. Реализует
`kodesChemistryModel` — химическую модель поверх `StandardChemistryModel`,
которая использует интегратор `Seulex` из библиотеки [KODES](src/KODES)
(CUDA) и кинетику, сгенерированную [pyJac](src/pyJacSource) (тоже CUDA).

Без CUDA библиотека тоже собирается (`wmake libso`) — просто химия считается
на CPU силами `StandardChemistryModel`, а GPU-часть (`src/KODES`,
`src/pyJacSource`) в сборку не включается.

## Требования

- OpenFOAM (испытано на v2412).
- CUDA Toolkit (`nvcc`, `nvlink`, `libcudart`/`libcusparse`/`libcublas`).
- `CUDA_DIR` в окружении должен указывать на корень CUDA Toolkit
  (`$CUDA_DIR/bin/nvcc`, `$CUDA_DIR/include`, `$CUDA_DIR/lib64`).
- **Компилятор OpenFOAM (`WM_COMPILER`) обязан быть `Nvidia`** (NVIDIA HPC
  SDK, `nvc++`), если библиотека собирается с `-cu`. Причина — см. ниже,
  раздел «Почему не GCC».

## Сборка

```bash
cd src
wmake libso            # без CUDA — химия считается на CPU
```

```bash
./Allwmake -cu          # с CUDA, архитектура по умолчанию (см. NVARCH=120 в Allwmake)
./Allwmake -cu 90       # с CUDA, явно указать compute capability (например, cc90)
./Allwmake -cl           # полная очистка (wclean + lnInclude + Make/files)
```

После смены `WM_COMPILER`, правки `wmake/*`-скриптов или переноса на другой
сервер **обязательно** делать чистую пересборку:

```bash
./Allwmake -cl
./Allwmake -cu
```
wmake не отслеживает изменения файлов `wmake/cuda` / `wmake/c++` как повод
пересобрать уже готовые `.o` (их пересборка триггерится только по `Make/options`
и исходникам), поэтому старые объектники могут молча остаться собранными
по старым флагам.

## Почему не GCC

`kodesChemistryModel` — шаблонный класс, который **инстанцируется прямо
внутри** `makeKodesChemistryModelTypes.C` и `makeKodesChemistrySolverTypes.C`
(макросы `makeKodesChemistryModelTypes`/`makeKodesChemistrySolverTypes`).
Эти `.C`-файлы транзитивно тянут `kodesChemistryModel.H` →
`KODES/include/Integrators/Seulex.cuh` → `KODES/src/Integrators/Seulex.cu`,
а там — реальный запуск CUDA-ядра:

```cpp
seulex_solve<ODESystem><<<this->blocks, this->threads, this->sharedMemSize>>>(...);
```

Синтаксис `<<<...>>>` понимают только `nvcc` и `nvc++ -cuda`. Обычный `g++`
падает с ошибкой парсинга — и это фундаментальное ограничение, не лечится
никакими `-I` или макросами. Поэтому:

- `WM_COMPILER=Nvidia` (`nvc++`) — работает, если добавить `-cuda -gpu=ccXX`
  (см. `wmake/c++`).
- `WM_COMPILER=Gcc`/`Clang`/etc. — работать не может без переработки кода
  (вынесение всех `__global__`/`__constant__`/`<<<>>>` конструкций за
  непрозрачный указатель, как это сделано в `AmgX4Foam` — там `.C`-файлы
  никогда не видят сырой CUDA-синтаксис, вся GPU-специфика спрятана в
  `.cu`-файлах за обычным C++ API).

`wmake/c++` сам решает, добавлять ли `-cuda -gpu=cc$(NVARCH)` к `CC`, по
значению `$(WM_COMPILER)` — руками ничего переключать не нужно, достаточно
поставить `WM_COMPILER=Nvidia` в окружении OpenFOAM перед сборкой:

```bash
export WM_COMPILER=Nvidia
source $WM_PROJECT_DIR/etc/bashrc
```

**Важно:** ту же переменную `WM_COMPILER=Nvidia` нужно засорсить и в сессии,
где потом запускается сам расчётный случай (solver) — иначе OpenFOAM будет
искать `libkodesChemistryModel.so` в `platforms/linux64GccDPInt32Opt/lib`,
а библиотека соберётся в `platforms/linux64NvidiaDPInt32Opt/lib`, и
`dlopen` не найдёт файл.

## Устройство CUDA-сборки

- `Allwmake -cu` собирает список исходников как `cat Make/filesCuda
  Make/filesIn > Make/files` и вызывает `wmake libso`.
- `wmake/cuda` подключается из `wmake/c++` только когда `have_cuda=true`,
  добавляет суффикс `.cu`, правило компиляции через `nvcc` и флаг
  `-rdc=true` (см. ниже, зачем).
- `wmake/link-cuda` **умышленно не подключён** (закомментирован в конце
  `wmake/cuda`): финальную линковку `.so` делает `nvc++` (через `LINKLIBSO`,
  собранный в `wmake/c++`), а не `nvcc`. Для связки «часть объектников из
  `nvcc -rdc=true`, часть — из `nvc++ -cuda`» этого достаточно: `nvc++`
  сам выполняет нужный device-link при финальной линковке.

### Пути к `.cu`-файлам в `Make/filesCuda`

wmake ищет файлы, перечисленные в `Make/files`, относительно `src/`. Все
исходники KODES и pyJac лежат не в `src/`, а в подпапках — пути в
`Make/filesCuda` должны быть полными относительными путями, а не голыми
именами файлов:

```
KODES/src/basic_linalg.cu
KODES/src/Resources/SeulexDeviceResources.cu
...
pyJacSource/grimech/out/chem_utils.cu
pyJacSource/grimech/out/jacobs/jacob_0.cu
pyJacSource/grimech/out/rates/rxn_rates_0.cu
...
```

При добавлении нового `.cu`-файла в KODES/pyJac — не забыть дописать его
путь (не имя!) в `src/Make/filesCuda`.

### Инклуды

`src/Make/options` явно перечисляет все директории с заголовками KODES и
pyJac (`-IKODES/include`, `-IKODES/include/Resources`, `-IpyJacSource/grimech/out`,
`-IpyJacSource/grimech/out/jacobs` и т.д.), плюс остаётся `-IlnInclude`
(плоский каталог симлинков, который строит `wmakeLnInclude`). Внутри самих
`.cu`/`.cuh`-файлов инклуды пишутся голыми именами (`#include
"basic_linalg.cuh"`) — резолвятся либо через явные `-I`, либо через
`lnInclude`. Дублей имён файлов в дереве нет — это важно, иначе плоский
`lnInclude` будет путать заголовки.

### `-rdc=true` и линковка device-кода

Сгенерированный pyJac-код активно вызывает `__device__`-функции между
разными `.cu`-файлами (`pyJacSystem.cu` → `dydt()` из `dydt.cu`, `jacob.cu`
→ `eval_jacob_0..10` из `jacobs/jacob_N.cu`, `rxn_rates.cu` → функции из
`rates/rxn_rates_N.cu`). По умолчанию `nvcc` компилирует каждый `.cu` в
изолированном режиме ("whole program") и не резолвит такие вызовы — падает
с `ptxas fatal: Unresolved extern function`. Отсюда `-rdc=true`
(relocatable device code) в `cuFLAGS` (`wmake/cuda`) — обязателен.

### `__constant__`/глобальные переменные — только `static`

Любая `__constant__`/`__device__` переменная, объявленная в заголовке
(`.cuh`), который подключается в **несколько** единиц трансляции, обязана
быть `static` (или C++17 `inline`). Иначе после включения `-rdc=true`
`nvlink` начинает по-настоящему сливать device-символы по всей программе и
падает с `Multiple definition of ...`. Наступили на эти грабли в
`KODES/include/Integrators/Seulex.cuh` (`absTol_`, `stepFactor*_`,
`nSeq_`, `cpu_`, `coeff_` и т.д. — все помечены `static`). При добавлении
новых `__constant__`-таблиц в заголовки — сразу ставить `static`.

### Финальный перевод строки

`wmkdepend` (сканер зависимостей OpenFOAM) падает с `parse error ...
perhaps missing a final newline`, если `.cu`-файл не заканчивается
переводом строки. Не фатально для сборки конкретного файла, но ломает
отслеживание зависимостей (пересборка по изменению заголовка перестаёт
работать). Проверять при коммите новых `.cu`-файлов.

## Диагностика типичных ошибок сборки

| Симптом | Причина | Что делать |
|---|---|---|
| `Нет правила для сборки цели .../*.cu.dep` | В `Make/filesCuda` указано голое имя файла, а не путь относительно `src/` | Прописать полный путь в `Make/filesCuda` |
| `wmkdepend: parse error ... missing a final newline` | В `.cu`-файле нет `\n` в конце | Добавить перевод строки в конец файла |
| `ptxas fatal: Unresolved extern function '...'` | Вызов `__device__`-функции из другого `.cu`-файла без `-rdc=true` | Проверить, что `-rdc=true` есть в `cuFLAGS` (`wmake/cuda`) |
| `nvlink error: Multiple definition of '...'` | `__constant__`/глобальная переменная объявлена (не `static`) в заголовке, подключённом в несколько `.o` | Пометить переменную `static` |
| `undefined symbol: __fatbinwrap_...` при `dlopen` рантайм-кейса | `.C`-файл, инстанцирующий CUDA-шаблоны (`kodesChemistryModel`), собран без `-cuda` (например, `nvc++` без флага, или `WM_COMPILER=Gcc`) | `WM_COMPILER=Nvidia` + `CC += -gpu=cc$(NVARCH) -cuda` (уже включено в `wmake/c++` автоматически при `WM_COMPILER=Nvidia`) |
| `fatal error: cuda/cmath: Нет такого файла` (или любой другой `#include <cuda_runtime.h>`/`__global__`/`<<<>>>` не парсится) | `.C`-файл собирается `g++` (`WM_COMPILER=Gcc`), который в принципе не понимает CUDA-синтаксис | Пересобрать с `WM_COMPILER=Nvidia`; на GCC текущая архитектура (сырой CUDA-код в заголовках, инклюдящихся в `.C`) не соберётся в принципе |
| `Could not load "libkodesChemistryModel.so"` при запуске кейса | Библиотека собрана под другим `WM_OPTIONS` (`WM_COMPILER`/precision/label-size), чем окружение солвера | Убедиться, что `WM_COMPILER` (и остальные `WM_*`) одинаковы при сборке и при запуске случая |
