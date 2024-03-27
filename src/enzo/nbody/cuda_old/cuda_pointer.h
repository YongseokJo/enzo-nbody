#include <assert.h>

template <typename T>
struct cudaPointer{
  T *dev_pointer;
  T *host_pointer;
  int size;

  cudaPointer(){
    dev_pointer  = NULL;
    host_pointer = NULL;
    size = 0;
  }
  //  ~cudaPointer(){
  // free();
  //  }
  void allocate(int _size){
		size = _size;
		void *p;
		CUDA_SAFE_CALL(cudaMalloc(&p, size * sizeof(T)));
		assert(p);
		dev_pointer = (T*)p;
		CUDA_SAFE_CALL(cudaMallocHost(&p, size * sizeof(T)));
		assert(p);
		host_pointer = (T*)p;
		/*
			 CUDA_SAFE_CALL(cudaMalloc(&dev_pointer, size * sizeof(T)));
			 CUDA_SAFE_CALL(cudaMallocHost(&host_pointer, size * sizeof(T)));
			 */
		assert(dev_pointer);
		assert(host_pointer);
  }

  void free(){
    CUDA_SAFE_CALL(cudaFree(dev_pointer));
    CUDA_SAFE_CALL(cudaFreeHost(host_pointer));
    dev_pointer  = NULL;
    host_pointer = NULL;
    size = 0;
  }

  void toDevice(int count){
    CUDA_SAFE_CALL(cudaMemcpy(dev_pointer, host_pointer, count * sizeof(T), cudaMemcpyHostToDevice));
  }

  void toDevice(){
    this->toDevice(size);
  }

  void toHost(int count){
    CUDA_SAFE_CALL(cudaMemcpy(host_pointer, dev_pointer, count * sizeof(T), cudaMemcpyDeviceToHost));
  }

  void toHost(){
    this->toHost(size);
  }

  T &operator [] (int i){
    return host_pointer[i];
  }

  operator T* (){
    return dev_pointer;
  }
};
