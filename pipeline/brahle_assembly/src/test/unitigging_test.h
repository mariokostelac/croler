#ifndef TEST_UNITIGGING_TEST_
#define TEST_UNITIGGING_TEST_

#include <overlap/read.h>
#include <overlap/overlap.h>
#include <layout/unitigging.h>

#include <vector>

namespace test {

class UnitiggingTest {
 public:
  UnitiggingTest();
  virtual ~UnitiggingTest();
  virtual bool run()=0;
  overlap::Read* makeRead(const char* data);
};

class UnitiggingIsTransitiveTest : public UnitiggingTest {
 public:
  UnitiggingIsTransitiveTest();
  virtual ~UnitiggingIsTransitiveTest();
  virtual bool run();
};

class UnitiggingContainmentTest : public UnitiggingTest {
 public:
  UnitiggingContainmentTest();
  virtual ~UnitiggingContainmentTest();
  virtual bool run();
};

class UnitiggingTransitiveTest : public UnitiggingTest {
 public:
  UnitiggingTransitiveTest();
  virtual ~UnitiggingTransitiveTest();
  virtual bool run();
};

class UnitiggingContigTest : public UnitiggingTest {
 public:
  UnitiggingContigTest();
  virtual ~UnitiggingContigTest();
  virtual bool run();
};

class UnitiggingCompleteTest : public UnitiggingTest {
  typedef overlap::Read* ReadPtr;
 public:
  explicit UnitiggingCompleteTest(
      const char* read_file,
      const char *overlap_file);
  virtual ~UnitiggingCompleteTest();
  virtual bool run();
 private:
  FILE *rfd_;
  FILE *ofd_;
};

class UnitiggingTestRunner {
 public:
  UnitiggingTestRunner();
  virtual ~UnitiggingTestRunner();

  void run();
  void addTest(UnitiggingTest* test);

 private:
  std::vector< UnitiggingTest* > tests_;
};

};  // namespace test

#endif  // TEST_UNITIGGING_TEST_
