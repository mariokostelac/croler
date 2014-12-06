#ifndef HASHLIST
#define HASHLIST
#include <cstdlib>
#include <vector>
#include <iterator>
#include <unordered_set>

template <typename K, typename T>
class HashList {
 public:
    class List {
     public:
         class iterator : std::iterator<std::input_iterator_tag, T> {
          public:
              explicit iterator(const HashList<K, T>& hashlist, K key, int pos)
                  : _hashlist(hashlist), _key(key), _pos(pos) {}

              iterator& operator++() {
                  do {
                      _pos = _hashlist._prev[_pos];
                  } while (_pos >= 0 && _hashlist._key[_pos] != _key);
                  return *this;
              }

              bool operator==(const iterator& rhs) {
                  return _pos == rhs._pos;
              }

              bool operator!=(const iterator& rhs) {
                  return _pos != rhs._pos;
              }

              const T& operator*() {
                  return _hashlist._value[_pos];
              }

          private:
              const HashList<K, T>& _hashlist;
              K _key;
              int _pos;
         };

         explicit List(const HashList& hl, K key) : _hashlist(hl), _key(key) {}

         K key() {
            return _key;
         }

         iterator begin() {
             int bucket = _key % _hashlist._bucket_num;
             int head = _hashlist._last[bucket];
             while (head >= 0 && _hashlist._key[head] != _key) {
                 head = _hashlist._prev[head];
             }
             return iterator(_hashlist, _key, head);
         }

         iterator end() {
             return iterator(_hashlist, _key, -1);
         }

     private:
         const HashList<K, T>& _hashlist;
         K _key;
    };

    explicit HashList(int bucket_num = 20000000) : _bucket_num(bucket_num), _size(0) {
        _last.resize(_bucket_num, -1);
    }

    void add(K key, T value) {
        int bucket = key % _bucket_num;
        _key.push_back(key);
        _value.push_back(value);
        _prev.push_back(_last[bucket]);
        _last[bucket] = _size;
        _keys.insert(key);
        ++_size;
    }

    const std::unordered_set<K>& keys() const {
        return _keys;
    }

    const List get_list(K key) const {
        return List(*this, key);
    }

    void shrink() {
        _value.shrink_to_fit();
        _prev.shrink_to_fit();
        _key.shrink_to_fit();
    }

    unsigned int size() {
        return _size;
    }

 private:
    unsigned int _bucket_num;
    unsigned int _size;
    std::vector<K> _key;
    std::vector<T> _value;
    std::vector<int> _prev;
    std::vector<int> _last;
    std::unordered_set<K> _keys;
};
#endif
