/*
 * Copyright (c) 2007 John Weaver
 * Copyright (c) 2015 Miroslav Krajicek
 *
 * This program is free software; you can redistribute it and/or modify
 * ... (Standard GPL Header) ...
 */

#include "lap.h"
#include <algorithm> // Needed for std::max

template<typename Data>
void LAPJV<Data>::solve(Matrix<Data> &m) {

  int original_rows = m.rows();
  int original_cols = m.columns();
  // LAPJV requires a square matrix. Pad to the largest dimension.
  int dim = std::max(original_rows, original_cols);

  // Helper lambda to seamlessly provide 0-cost padding for out-of-bounds access
  auto cost_at = [&](int r, int c) -> Data {
      if (r < original_rows && c < original_cols) {
        return m(r, c);
      }
      return 0; // 0 cost ensures "dummy" assignments don't penalize the real assignments
  };

  bool unassignedfound;
  row i, imin, numfree = 0, prvnumfree, f, i0, k, freerow, *pred, *freeunassigned;
  col j, j1, j2, endofpath, last, low, up, *collist, *matches;
  cost min, h, umin, usubmin, v2, *d;

  col *rowsol = new col[dim];
  row *colsol = new row[dim];
  cost *u = new cost[dim];
  cost *v = new cost[dim];

  freeunassigned = new row[dim];
  collist = new col[dim];
  matches = new col[dim];
  d = new cost[dim];
  pred = new row[dim];

  for (i = 0; i < dim; i++) matches[i] = 0;

  // COLUMN REDUCTION
  for (j = dim; j--;)
  {
    min = cost_at(0, j);
    imin = 0;
    for (i = 1; i < dim; i++)
      if (cost_at(i, j) < min) {
        min = cost_at(i, j);
        imin = i;
      }
    v[j] = min;
    if (++matches[imin] == 1) {
      rowsol[imin] = j;
      colsol[j] = imin;
    } else if (v[j] < v[rowsol[imin]]) {
      int j1 = rowsol[imin];
      rowsol[imin] = j;
      colsol[j] = imin;
      colsol[j1] = -1;
    } else
      colsol[j] = -1;
  }

  // REDUCTION TRANSFER
  for (i = 0; i < dim; i++)
    if (matches[i] == 0)
      freeunassigned[numfree++] = i;
    else if (matches[i] == 1)
    {
      j1 = rowsol[i];
      min = BIG;
      for (j = 0; j < dim; j++)
        if (j != j1)
          if (cost_at(i, j) - v[j] < min) min = cost_at(i, j) - v[j];
      v[j1] = v[j1] - min;
    }

  // AUGMENTING ROW REDUCTION
  int loopcnt = 0;
  do {
    loopcnt++;

    k = 0;
    prvnumfree = numfree;
    numfree = 0;
    while (k < prvnumfree) {
      i = freeunassigned[k];
      k++;

      umin = cost_at(i, 0) - v[0];
      j1 = 0;
      usubmin = BIG;
      for (j = 1; j < dim; j++) {
        h = cost_at(i, j) - v[j];
        if (h < usubmin)
          if (h >= umin) {
            usubmin = h;
            j2 = j;
          } else {
            usubmin = umin;
            umin = h;
            j2 = j1;
            j1 = j;
          }
      }

      i0 = colsol[j1];
      if (umin < usubmin)
        v[j1] = v[j1] - (usubmin - umin);
      else
      if (i0 > -1)
      {
        j1 = j2;
        i0 = colsol[j2];
      }

      rowsol[i] = j1;
      colsol[j1] = i;

      if (i0 > -1)
        if (umin < usubmin)
          freeunassigned[--k] = i0;
        else
          freeunassigned[numfree++] = i0;
    }
  } while (loopcnt < 2);

  // AUGMENT SOLUTION for each free row.
  for (f = 0; f < numfree; f++) {
    freerow = freeunassigned[f];

    for (j = dim; j--;) {
      d[j] = cost_at(freerow, j) - v[j];
      pred[j] = freerow;
      collist[j] = j;
    }

    low = 0;
    up = 0;
    unassignedfound = false;
    do {
      if (up == low)
      {
        last = low - 1;

        min = d[collist[up++]];
        for (k = up; k < dim; k++) {
          j = collist[k];
          h = d[j];
          if (h <= min) {
            if (h < min)
            {
              up = low;
              min = h;
            }
            collist[k] = collist[up];
            collist[up++] = j;
          }
        }
        for (k = low; k < up; k++)
          if (colsol[collist[k]] < 0) {
            endofpath = collist[k];
            unassignedfound = true;
            break;
          }
      }

      if (!unassignedfound) {
        j1 = collist[low];
        low++;
        i = colsol[j1];
        h = cost_at(i, j1) - v[j1] - min;

        for (k = up; k < dim; k++) {
          j = collist[k];
          v2 = cost_at(i, j) - v[j] - h;
          if (v2 < d[j]) {
            pred[j] = i;
            if (v2 == min)
              if (colsol[j] < 0) {
                endofpath = j;
                unassignedfound = true;
                break;
              }
              else {
                collist[k] = collist[up];
                collist[up++] = j;
              }
            d[j] = v2;
          }
        }
      }
    } while (!unassignedfound);

    for (k = last + 1; k--;) {
      j1 = collist[k];
      v[j1] = v[j1] + d[j1] - min;
    }

    do {
      i = pred[endofpath];
      colsol[endofpath] = i;
      j1 = endofpath;
      endofpath = rowsol[i];
      rowsol[i] = j1;
    } while (i != freerow);
  }

  // calculate optimal cost.
  cost lapcost = 0;
  for (i = dim; i--;) {
    j = rowsol[i];
    u[i] = cost_at(i, j) - v[j];
    lapcost = lapcost + cost_at(i, j);
  }

  // Write assignments back into the input matrix
  for (size_t r = 0; r < m.rows(); ++r) {
    for (size_t c = 0; c < m.columns(); ++c) {
      if (rowsol[r] == (col)c) {
        m(r, c) = 0.0;
      } else {
        m(r, c) = -1.0;
      }
    }
  }

  // free reserved memory.
  delete[] pred;
  delete[] freeunassigned;
  delete[] collist;
  delete[] matches;
  delete[] d;

  delete[] rowsol;
  delete[] colsol;
  delete[] u;
  delete[] v;
}
template class LAPJV<double>;
template class LAPJV<float>;
template class LAPJV<int>;
