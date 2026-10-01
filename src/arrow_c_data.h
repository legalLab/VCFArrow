// src/arrow_c_data.h
//
// Arrow C Data Interface structs (ABI-stable, copied from the Arrow spec:
// https://arrow.apache.org/docs/format/CDataInterface.html) and helpers for
// passing them between R (the arrow package) and C++ as external pointers.

#ifndef VCFARROW_ARROW_C_DATA_H
#define VCFARROW_ARROW_C_DATA_H

#include <Rinternals.h>
#include <cstdint>

#ifndef ARROW_C_DATA_INTERFACE
#define ARROW_C_DATA_INTERFACE

#define ARROW_FLAG_DICTIONARY_ORDERED 1
#define ARROW_FLAG_NULLABLE 2
#define ARROW_FLAG_MAP_KEYS_SORTED 4

extern "C" {
struct ArrowSchema {
  const char* format;
  const char* name;
  const char* metadata;
  int64_t flags;
  int64_t n_children;
  struct ArrowSchema** children;
  struct ArrowSchema* dictionary;
  void (*release)(struct ArrowSchema*);
  void* private_data;
};

struct ArrowArray {
  int64_t length;
  int64_t null_count;
  int64_t offset;
  int64_t n_buffers;
  int64_t n_children;
  const void** buffers;
  struct ArrowArray** children;
  struct ArrowArray* dictionary;
  void (*release)(struct ArrowArray*);
  void* private_data;
};
}

#endif  // ARROW_C_DATA_INTERFACE

// Finalizers for external pointers owning an ArrowArray / ArrowSchema struct:
// release its contents if nobody took them over, then free the struct.
inline void arrow_array_finalizer(SEXP x) {
  ArrowArray* a = static_cast<ArrowArray*>(R_ExternalPtrAddr(x));
  if (a) {
    if (a->release) a->release(a);
    delete a;
    R_ClearExternalPtr(x);
  }
}

inline void arrow_schema_finalizer(SEXP x) {
  ArrowSchema* s = static_cast<ArrowSchema*>(R_ExternalPtrAddr(x));
  if (s) {
    if (s->release) s->release(s);
    delete s;
    R_ClearExternalPtr(x);
  }
}

#endif  // VCFARROW_ARROW_C_DATA_H
