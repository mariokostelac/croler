#include <queue>
#include <list>

template <typename T>
class FixedMinQueue {
 public:
     explicit FixedMinQueue(int n) {
         size = n;
         q = new std::list<T>();
         elems = new std::queue<T>();
     }

     ~FixedMinQueue() {
         if (q != 0)        delete(q);
         if (elems != 0)    delete(elems);
     }

     void push(T elem) {
         while (!q->empty() && q->back() > elem) q->pop_back();
         q->push_back(elem);
         elems->push(elem);
         if (elems->size() > size) pop();
     }

     T pop() {
         T ret = elems->front();
         if (q->front() == ret) q->pop_front();
         elems->pop();
         return ret;
     }

     T min() {
         return q->front();
     }

     void clear() {
         q->clear();
         while (elems->size()) elems->pop();
     }

 private:
     int size;
     std::queue<T>* elems;
     std::list<T>* q;
};

