use crate::state::*;

impl Engine {
    pub(crate) fn heap_replace_top(&mut self, key: f64, value: usize) {
        let num = self.heap.len();
        let mut i = 0;
        let mut j1 = 1;
        let mut j2 = 2;
        if num > 2 {
            loop {
                let j = if self.heap[j2].key <= self.heap[j1].key {
                    j2
                } else {
                    j1
                };
                if key <= self.heap[j].key {
                    j1 = num;
                    j2 = num;
                } else {
                    self.heap[i] = self.heap[j];
                    j1 = 2 * j + 1;
                    j2 = 2 * j + 2;
                    i = j;
                }
                if j2 >= num {
                    break;
                }
            }
        }
        if j1 < num && self.heap[j1].key < key {
            self.heap[i] = self.heap[j1];
            self.heap[j1] = HeapEntry { key, value };
            return;
        }
        self.heap[i] = HeapEntry { key, value };
    }

    pub(crate) fn heap_insert(&mut self, key: f64, value: usize) {
        let mut j = self.heap.len();
        self.heap.push(HeapEntry { key, value });
        while j > 0 {
            let parent = (j - 1) / 2;
            if self.heap[parent].key <= key {
                break;
            }
            self.heap[j] = self.heap[parent];
            j = parent;
        }
        self.heap[j] = HeapEntry { key, value };
    }

    pub(crate) fn heap_push(&mut self, key: f64, value: usize) {
        if self.heap_top_pending {
            self.heap_replace_top(key, value);
        } else {
            self.heap_insert(key, value);
        }
        self.heap_top_pending = false;
        self.heap_top_valid = false;
    }

    /// Lazy pop: the returned root stays in place until the next push or pop replaces it.
    pub(crate) fn heap_pop_min(&mut self) -> Option<(f64, usize)> {
        if self.heap_top_pending && self.heap_top_valid {
            if self.heap.len() <= 1 {
                self.heap.clear();
                return None;
            }
            let last = self.heap.pop().expect("heap holds more than one entry");
            self.heap_replace_top(last.key, last.value);
            self.heap_top_valid = false;
            self.heap_top_pending = false;
        } else if self.heap.is_empty() {
            return None;
        }
        let top = self.heap[0];
        self.heap_top_valid = true;
        self.heap_top_pending = true;
        Some((top.key, top.value))
    }
}
