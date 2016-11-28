/**
 @file Find_k_max.h
 Find the k greatest elements out of n.
 */

#ifndef FIND_K_MAX_H
#define FIND_K_MAX_H


/**
 * \brief Swap parameters.
 *
 * @param a - first parameter
 * @param b - second parameter
 */
void swap(float &a, float &b) {
  a = a + b;
  b = a - b;
  a = a - b;
}

/**
 * \brief Swap parameters.
 *
 * @param a - first parameter
 * @param b - second parameter
 */
void swap(size_t& a, size_t& b) {
  a = a + b;
  b = a - b;
  a = a - b;
}

/**
 * \brief Make min heap.
 *
 * @param a         - the data
 * @param size      - size of data
 * @param i         - helping index
 * @param indices   - indices of the data
 */
void minHeapify(float a[], int size, int i, size_t indices[]) {
  int l = 2 * i;
  int r = 2 * i + 1;
  int smallest = i;
  if (l < size && a[l] < a[smallest])
    smallest = l;
  if (r < size && a[r] < a[smallest])
    smallest = r;
  if (smallest != i) {
    swap(a[i], a[smallest]);
    swap(indices[i], indices[smallest]);
    minHeapify(a, size, smallest, indices);
  }

}

/**
 * \brief Build min heap.
 *
 * @param a         - the data
 * @param size      - size of data
 * @param indices   - indices of the data
 */
void buildMinHeap(float a[], int size, size_t indices[]) {
  for (int i = size / 2; i >= 0; i--)
    minHeapify(a, size, i, indices);
}

/**
 * \brief Find the k max elements.
 *
 * @param a         - the data
 * @param size      - size of data
 * @param k         - number of max elements
 * @param indices   - indices of the data
 */
void kthLargest(float a[], int size, int k, size_t indices[]) {
  float minHeap[k];
  int i;
  for (i = 0; i < k; i++) {
    minHeap[i] = a[i];
    indices[i] = i;
  }
  buildMinHeap(minHeap, k, indices);
  for (i = k; i < size; i++) {
    if (a[i] > minHeap[0]) {
      minHeap[0] = a[i];
      indices[0] = i;
      minHeapify(minHeap, k, 0, indices);
    }
  }
}

#endif
