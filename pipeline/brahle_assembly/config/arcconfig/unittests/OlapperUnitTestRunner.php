<?php

class OlapperUnitTestRunner extends ArcanistBaseUnitTestEngine {
  public function run() {
    $future = new ExecFuture("make test");
    $future->setCWD($this->getWorkingCopy()->getProjectRoot());
    list($stdout, $stderr) = $future->resolvex();
    return $this->parseTestResults($stderr);
  }

  protected function parseTestResults($data) {
    $tests = explode("\n", $data);
    $results = array();
    foreach ($tests as $test) {
      $result = $this->parseOneTestResult($test);
      if ($result !== null) {
        $results[] = $result;
      }
    }
    return $results;
  }

  protected function parseOneTestResult($data) {
    $data = $this->removeCLIColors($data);
    $result = new ArcanistUnitTestResult();
    $pattern = "/";
    $pattern .= "(?P<test_result>\w+):";
    $pattern .= " (?P<name>\w+)";
    $pattern .= " \((?P<time>\d+.\d+s)\)";
    $pattern .= "/";
    if (!preg_match($pattern, $data, $matches)) {
      return null;
    }

    if ($matches["test_result"] == "OK") {
      $test_result = ArcanistUnitTestResult::RESULT_PASS;
    } else {
      $test_result = ArcanistUnitTestResult::RESULT_FAIL;
    }
    $name = $this->demangleName($matches["name"]);
    $duration = floatval($matches["time"]);

    $result->setName($name);
    $result->setResult($test_result);
    $result->setDuration($duration);
    return $result;
  }

  protected function removeCLIColors($data) {
    $data = preg_replace("/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]/", "", $data);
    return $data;
  }

  protected function demangleName($name) {
    //TODO: Implement this
    return $name;
  }
}


