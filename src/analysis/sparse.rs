//! Defines matrices in coordinate format (COO) and Compressed Sparse Row (CSR).

use crate::error::FemError;

/// A sparse matrix in Coordinate List format.
///
/// COO stores matrix entries in three parallel vectors:
///
/// rows[k] - row index of entry k
/// cols[k] - column index of entry k
/// values[k] - value of entry k
#[derive(Debug, Clone, PartialEq)]
pub struct CooMatrix {
    nrows: usize,
    ncols: usize,
    rows: Vec<usize>,
    cols: Vec<usize>,
    values: Vec<f64>,
}

impl CooMatrix {
    /// Creates an empty COO matrix with the requested dimensions.
    pub fn new(nrows: usize, ncols: usize) -> Self {
        Self { nrows, ncols, rows: Vec::new(), cols: Vec::new(), values: Vec::new() }
    }

    /// Returns the number of matrix rows.
    pub fn nrows(&self) -> usize {
        self.nrows
    }

    /// Returns the number of matrix columns.
    pub fn ncols(&self) -> usize {
        self.ncols
    }

    /// Returns the number of stored entries.
    pub fn nnz(&self) -> usize {
        self.values.len()
    }

    /// Adds one matrix entry to the COO representation.
    pub fn push(&mut self, row: usize, col: usize, value: f64) -> Result<(), FemError> {
        if row >= self.nrows {
            return Err(FemError::InvalidDofIndex { index: row });
        }

        if col >= self.ncols {
            return Err(FemError::InvalidDofIndex { index: col });
        }

        if value == 0.0 {
            return Ok(());
        }

        self.rows.push(row);
        self.cols.push(col);
        self.values.push(value);

        Ok(())
    }

    /// Sums entries that have the same `(row, column)` position.
    pub fn sum_duplicates(&mut self) {
        let rows = std::mem::take(&mut self.rows);
        let cols = std::mem::take(&mut self.cols);
        let values = std::mem::take(&mut self.values);

        let mut entries: Vec<(usize, usize, f64)> =
            rows.into_iter().zip(cols).zip(values).map(|((row, col), value)| (row, col, value)).collect();

        entries.sort_unstable_by(|left, right| left.0.cmp(&right.0).then(left.1.cmp(&right.1)));

        let mut entries = entries.into_iter();

        let Some((mut current_row, mut current_col, mut current_value)) = entries.next() else {
            return;
        };

        for (row, col, value) in entries {
            if row == current_row && col == current_col {
                current_value += value;
                continue;
            }

            if current_value != 0.0 {
                self.rows.push(current_row);
                self.cols.push(current_col);
                self.values.push(current_value);
            }

            current_row = row;
            current_col = col;
            current_value = value;
        }

        if current_value != 0.0 {
            self.rows.push(current_row);
            self.cols.push(current_col);
            self.values.push(current_value);
        }
    }

    /// Converts the COO matrix into CSR format.
    pub fn into_csr(mut self) -> CsrMatrix {
        self.sum_duplicates();

        let mut row_offsets = vec![0; self.nrows + 1];

        // Count how many entries belong to every row.
        for &row in &self.rows {
            row_offsets[row + 1] += 1;
        }

        // Convert row counts into cumulative offsets.
        for row in 0..self.nrows {
            row_offsets[row + 1] += row_offsets[row];
        }

        CsrMatrix { nrows: self.nrows, ncols: self.ncols, row_offsets, column_indices: self.cols, values: self.values }
    }
}

/// A sparse matrix in Compressed Sparse Row format.
///
/// CSR uses three arrays:
///
/// - `row_offsets` has length `nrows + 1`;
/// - `column_indices` stores the column of every non-zero entry;
/// - `values` stores the corresponding matrix values.
#[derive(Debug, Clone, PartialEq)]
pub struct CsrMatrix {
    nrows: usize,
    ncols: usize,
    row_offsets: Vec<usize>,
    column_indices: Vec<usize>,
    values: Vec<f64>,
}

impl CsrMatrix {
    /// Returns the number of matrix rows.
    pub fn nrows(&self) -> usize {
        self.nrows
    }

    /// Returns the number of matrix columns.
    pub fn ncols(&self) -> usize {
        self.ncols
    }

    /// Returns the number of stored non-zero entries.
    pub fn nnz(&self) -> usize {
        self.values.len()
    }

    /// Computes `y = A * x`.
    pub fn mul_vector(&self, x: &[f64], y: &mut [f64]) -> Result<(), FemError> {
        if x.len() != self.ncols {
            return Err(FemError::InvalidVectorLength { vector: "input", expected: self.ncols, actual: x.len() });
        }

        if y.len() != self.nrows {
            return Err(FemError::InvalidVectorLength { vector: "output", expected: self.nrows, actual: y.len() });
        }

        for (row, output) in y.iter_mut().enumerate().take(self.nrows) {
            let start = self.row_offsets[row];
            let end = self.row_offsets[row + 1];

            let mut row_sum = 0.0;

            for entry_index in start..end {
                let column = self.column_indices[entry_index];
                row_sum += self.values[entry_index] * x[column];
            }

            *output = row_sum;
        }

        Ok(())
    }

    /// Returns the value at `(row, column)`.
    pub fn value_at(&self, row: usize, column: usize) -> Result<f64, FemError> {
        if row >= self.nrows {
            return Err(FemError::InvalidDofIndex { index: row });
        }

        if column >= self.ncols {
            return Err(FemError::InvalidDofIndex { index: column });
        }

        let start = self.row_offsets[row];
        let end = self.row_offsets[row + 1];

        match self.column_indices[start..end].binary_search(&column) {
            Ok(relative_index) => Ok(self.values[start + relative_index]),
            Err(_) => Ok(0.0),
        }
    }

    /// Sets the value at `(row, column)`.
    pub fn set_value(&mut self, row: usize, column: usize, value: f64) -> Result<(), FemError> {
        if row >= self.nrows {
            return Err(FemError::InvalidDofIndex { index: row });
        }

        if column >= self.ncols {
            return Err(FemError::InvalidDofIndex { index: column });
        }

        let start = self.row_offsets[row];
        let end = self.row_offsets[row + 1];

        match self.column_indices[start..end].binary_search(&column) {
            Ok(relative_index) => {
                let index = start + relative_index;

                if value == 0.0 {
                    self.column_indices.remove(index);
                    self.values.remove(index);

                    for offset in (row + 1)..=self.nrows {
                        self.row_offsets[offset] -= 1;
                    }
                } else {
                    self.values[index] = value;
                }
            }
            Err(relative_index) => {
                if value == 0.0 {
                    return Ok(());
                }

                let index = start + relative_index;

                self.column_indices.insert(index, column);
                self.values.insert(index, value);

                for offset in (row + 1)..=self.nrows {
                    self.row_offsets[offset] += 1;
                }
            }
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::{CooMatrix, CsrMatrix};
    use crate::error::FemError;

    #[test]
    fn sums_duplicate_entries() {
        let mut matrix = CooMatrix::new(3, 3);

        matrix.push(1, 2, 2.0).unwrap();
        matrix.push(1, 2, 3.0).unwrap();
        matrix.push(0, 0, 4.0).unwrap();

        matrix.sum_duplicates();

        assert_eq!(matrix.nnz(), 2);
        assert_eq!(matrix.rows, vec![0, 1]);
        assert_eq!(matrix.cols, vec![0, 2]);
        assert_eq!(matrix.values, vec![4.0, 5.0]);
    }

    #[test]
    fn removes_entries_that_cancel_to_zero() {
        let mut matrix = CooMatrix::new(2, 2);

        matrix.push(0, 1, 7.0).unwrap();
        matrix.push(0, 1, -7.0).unwrap();

        matrix.sum_duplicates();

        assert_eq!(matrix.nnz(), 0);
    }

    #[test]
    fn converts_coo_to_csr_and_multiplies_vector() {
        let mut matrix = CooMatrix::new(3, 3);

        matrix.push(0, 0, 2.0).unwrap();
        matrix.push(0, 2, 1.0).unwrap();
        matrix.push(1, 1, 3.0).unwrap();
        matrix.push(2, 0, 4.0).unwrap();

        let csr = matrix.into_csr();

        let mut result = vec![0.0; 3];
        csr.mul_vector(&[1.0, 2.0, 3.0], &mut result).unwrap();

        assert_eq!(result, vec![5.0, 6.0, 4.0]);
    }

    #[test]
    fn validates_vector_dimensions() {
        let mut matrix = CooMatrix::new(2, 2);
        matrix.push(0, 0, 1.0).unwrap();

        let csr = matrix.into_csr();
        let mut result = vec![0.0; 2];

        let error = csr.mul_vector(&[1.0], &mut result);

        assert!(matches!(error, Err(FemError::InvalidVectorLength { vector: "input", expected: 2, actual: 1 })));
    }

    #[test]
    fn rejects_invalid_matrix_indices() {
        let mut matrix = CooMatrix::new(2, 2);

        let error = matrix.push(2, 0, 1.0);

        assert!(matches!(error, Err(FemError::InvalidDofIndex { index: 2 })));
    }

    #[test]
    fn supports_empty_rows_in_csr() {
        let mut matrix = CooMatrix::new(3, 3);
        matrix.push(1, 2, 5.0).unwrap();

        let csr: CsrMatrix = matrix.into_csr();

        assert_eq!(csr.row_offsets, vec![0, 0, 1, 1]);
    }
}
