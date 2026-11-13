# MATLAB cheat sheet – VSM167 / BasicFEM

Short reference for the homework and computer assignments.
`A`, `K` = matrices, `v`, `a`, `f` = vectors, `x` = scalar.

---

## 1. Paths (top of every script)

```matlab
here = fileparts(mfilename('fullpath'));              % folder of this script
if ~isfolder(fullfile(here, '..', '..', 'src'))       % Run Section uses a temp copy
    here = fileparts(matlab.desktop.editor.getActiveFilename);
end
addpath(here)                                         % your own functions
addpath(genpath(fullfile(here, '..', '..', 'src')))   % BasicFEM src + subfolders
```

| Rule | Why |
|---|---|
| Never hardcode `C:\Users\...` | Breaks on every other computer |
| Build file paths with `fullfile(here, ...)` | Correct separator on Windows/Mac/Linux |
| Function name = file name | MATLAB finds functions by **file** name |
| Code before `function` turns the file into a script | Keep function files starting with `function` |
| `addpath` only lasts for the MATLAB session | Each script must add its own paths |
| Don't use `savepath` for course work | The grader won't have your saved path |

Check: `which my_function` shows which file MATLAB will use.
Test like a grader: `restoredefaultpath`, `cd(userpath)`, then `run(...)` the script.

---

## 2. Matrices

### Creating

```matlab
A  = zeros(3, 3);            % 3x3 of zeros (preallocate before loops)
v  = zeros(size(xi));        % same shape as xi
I  = eye(3);                 % identity
x  = linspace(-1, 1, 101);   % 101 points from -1 to 1 (row vector)
k  = 1:5;                    % [1 2 3 4 5]
c  = [1, 0, 0]';             % ' transposes -> column vector (3x1)
A  = [c1, c2, c3];           % columns side by side  -> 3x3
A  = [r1; r2; r3];           % rows stacked (;)      -> 3x?
```

### Size

| Code | Gives |
|---|---|
| `size(A)` | `[rows cols]` |
| `size(A, 1)` | number of rows |
| `size(A, 2)` | number of columns |
| `[r, c] = size(A)` | both |
| `numel(A)` | total number of elements |
| `length(A)` | largest dimension |

### Selecting – `A(row, column)`, `:` = "all"

| Code | Selects |
|---|---|
| `A(2, 3)` | one element |
| `A(2, :)` | row 2 (1×n) |
| `A(:, 3)` | column 3 (n×1) |
| `A(:, end)` | last column |
| `A(:, [1 3])` | columns 1 and 3 |
| `A(1:2, 2:3)` | sub-block |
| `A(:)` | everything as one long column (column by column) |
| `v(v > 0)` | elements where the condition is true (logical indexing) |

> Two indexings in a row are **not** allowed: `A(:, j)(1)` ✘ → use `A(1, j)` or `aj = A(:, j); aj(1)` ✔

Loop over the columns of a matrix:

```matlab
for k = 1:size(A, 2)
    a = A(:, k);     % k-th column
end
```

### Multiplication

| Code | Type | Rule |
|---|---|---|
| `A * B` | matrix product | inner sizes must match: (m×**n**)(**n**×p) → m×p |
| `A .* B` | element by element | same size (or expandable) |
| `A .^ 2`, `A ./ B` | element by element | note the dot |
| `x * A` | scalar times matrix | always fine |
| `A'` | transpose | rows ↔ columns |

The order matters with vectors (`N` is a 1×3 row):

| Code | Size | Meaning |
|---|---|---|
| `N * a` (a = 3×1) | 1×1 | $\sum_i N_i a_i$ – value of the approximation |
| `N' * N` | 3×3 | all pairs $N_i N_j$ (outer product) → `K` |
| `N * N'` | 1×1 | $\sum_i N_i^2$ – usually **not** what you want |
| `N' * g` | 3×1 | $N_i\,g$ → `f` |

Many points at once: `Nmat` (101×3) times `A` (3×6) gives a 101×6 matrix, with one curve per column.

**Implicit expansion:** a row + a column (or a scalar + a matrix) silently expands to a full matrix.
`[1 2 3] + [1; 2]` gives a 2×3 matrix and **no error**. When results are oddly sized, check `size(...)`.

### Solving and checking

```matlab
a = K \ f;                    % solve K a = f   (don't use inv(K)*f)
rank(K)                       % must equal size(K,1), otherwise K is singular
norm(K - K')                  % ~0 if K is symmetric
abs(x - x_ref) < 1e-10        % compare floats with a tolerance, never with ==
```

The warning *"Matrix is singular to working precision"* or `NaN`s in the result means `K` isn't invertible (e.g. too few quadrature points).

### Sparse vs full

- Small element matrices (3×3, 8×8): use `zeros` (full).
- Large global stiffness matrices (thousands of DoFs, mostly zeros): use `sparse(ndofs, ndofs)`.
- Sparse + full = **full**, so a sparse 3×3 gains nothing.

---

## 3. Plotting

### Basic template

```matlab
figure                                   % new window, so nothing is drawn on top of older plots
hold on                                  % keep adding curves to the same axes
plot(x, y1, 'k--', 'LineWidth', 2)       % black dashed
plot(x, y2, 'LineWidth', 1.5)
hold off
grid on
xlabel('\xi'), ylabel('f(\xi; a)')       % \xi = Greek letter (TeX)
title('Quadratic line element')
legend('e^\xi', '3 points', '4 points', 'Location', 'best')
```

- `plot(x, Y)` with `Y` as a matrix draws **one curve per column**, so you don't need a loop or `hold on`.
- `x` and `y` must have the same number of points. `plot(x, scalar)` shows nothing useful.
- Write name/value options as `'LineWidth', 2`. The form `LineWidth=2` only works in R2021a+.

### Line styles

| Code | Result |
|---|---|
| `'-'` `'--'` `':'` `'-.'` | solid, dashed, dotted, dash-dot |
| `'o'` `'x'` `'s'` `'.'` | markers |
| `'r'` `'b'` `'k'` `'g'` | red, blue, black, green |
| `'ro-'` | red line with circle markers |

### Useful extras

```matlab
plot(x, y, 'DisplayName', 'a = [1 0 0]')   % label per curve ...
legend show                                 % ... then show all labels
xline(0, ':')                               % vertical reference line
yline(1, ':')                               % horizontal reference line
axis equal                                  % same scale on both axes (meshes!)
xlim([-1 1]), ylim([-0.2 1.2])              % axis limits

tiledlayout(1, 2)                           % several plots in one figure
nexttile, plot(x, y1), title('Fit')
nexttile, plot(x, y1 - y2), title('Error')
```

LaTeX in labels:
`xlabel('$\xi$', 'Interpreter', 'latex')`

### Saving figures (assignments → `resultfiles`)

```matlab
exportgraphics(gcf, fullfile(here, 'resultfiles', 'fit_exp.png'), 'Resolution', 300)
exportgraphics(gcf, fullfile(here, 'resultfiles', 'fit_exp.pdf'))   % vector, good for reports
```

---

## 4. Printing

| Code | Use |
|---|---|
| leave out the `;` | quick look: `K` prints the name + value |
| `disp(K)` | value only, no name |
| `fprintf(...)` | formatted output, for reports and checks |
| `sprintf(...)` | same, but returns a string (for titles/legends) |
| `mat2str(a', 4)` | vector/matrix as text, 4 significant digits, e.g. `'[0.5 1 0.75]'` |
| `format long` / `format short` | more or fewer digits on screen |

### `fprintf` formats

| Format | Prints |
|---|---|
| `%d` | integer |
| `%.4f` | 4 decimals: `2.6285` |
| `%10.4f` | width 10, 4 decimals (aligns columns) |
| `%.3e` | scientific: `1.234e-05` |
| `%s` | text |
| `\n` | new line (fprintf **doesn't** add one itself) |

```matlab
fprintf('nqp = %d: a = [%.4f, %.4f, %.4f]\n', nqp, a)
fprintf('max error = %.3e\n', max(abs(u - exp(xi))))
```

Print a matrix **row by row** (fprintf reads column by column, so transpose):

```matlab
fprintf('%10.4f %10.4f %10.4f\n', K')
```

Dynamic titles/legends:

```matlab
title(sprintf('Fit with %d quadrature points', nqp))
legend(sprintf('a = %s', mat2str(a', 3)))
```

Clear the screen: `clc`. Clear variables: `clear` (this **doesn't** clear the path).

---

## 5. Common errors seen so far

| Message | Usual cause |
|---|---|
| `Undefined function 'xyz' for input arguments of type 'double'` | Folder not on the path, wrong current folder, or the file name doesn't match the function name |
| `Unrecognized function or variable 'here'` | Variable used before it's defined (only worked because of an old workspace) |
| `Index exceeds the number of array elements` | Asking for `N(3)` from a 1×2 result, or a row/column mix-up (`a(j,:)` vs `a(:,j)`) |
| `Invalid array indexing` | Double indexing like `A(:,j)(1)` |
| `Unable to perform assignment because the left and right sides have a different number of elements` | Product in the wrong order, e.g. `a * N` (3×3) instead of `N * a` (1×1) |
| `Matrix is singular to working precision` / `NaN` | `K` not invertible (too few quadrature points, missing boundary conditions) |
| Only one curve shows | Missing `hold on` |
| New plot drawn on top of an old one | Missing `figure` |
| Flat line with one spike | Writing `u(i)` with the wrong loop counter |

**Habits:** give loop counters meaningful names (`nqp`, `k`, `iel`), preallocate with `zeros`, and check `size(...)` whenever something looks off.
