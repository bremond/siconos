(list (channel
        (name 'guix)
        (url "https://git.savannah.gnu.org/git/guix.git")
        (branch "master")
        (commit
          "1b2505217cf222d98cc960b8510660976a01cfa1")
        (introduction
          (make-channel-introduction
            "9edb3f66fd807b096b48283debdcddccfea34bad"
            (openpgp-fingerprint
              "BBB0 2DDF 2CEA F6A8 0D1D  E643 A2A0 6DF2 A33A 54FA"))))
      (channel
        (name 'guix-hpc)
        (url "https://gitlab.inria.fr/bremond/guix-hpc.git")
        (branch "add-siconos")
        (commit
          "05b87e068001a7b00aa3b61eb21c89f83335c116")))
