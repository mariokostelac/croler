<?php

class OlapperLintEngine extends ArcanistLintEngine {
  public function buildLinters() {
    $linters = array();
    $paths = $this->getPaths();

    foreach ($paths as $key => $path) {
      if (preg_match('@^lib/@', $path)) {
        unset($paths[$key]);
      }
      if (preg_match('@^bin/@', $path)) {
        unset($paths[$key]);
      }
    }

    $cpp_paths = preg_grep("/\.(cpp|cxx|cc|h|hpp)$/", $paths);
    $linters[] = id(new ArcanistCpplintLinter())->setPaths($cpp_paths);
    $linters[] = id(new ArcanistFilenameLinter())->setPaths($paths);

    return $linters;
  }

}

