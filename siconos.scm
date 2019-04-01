(define-module (siconos)
  #:use-module (guix packages)
  #:use-module (guix download)
  #:use-module (guix git-download)
  #:use-module (guix build-system cmake)
  #:use-module (guix licenses)
  #:use-module (gnu packages swig)
  #:use-module (gnu packages boost)
  #:use-module (gnu packages swig)
  #:use-module (gnu packages python)
  #:use-module (gnu packages python-xyz)
  #:use-module (gnu packages multiprecision)
  #:use-module (gnu packages gcc)
  #:use-module (gnu packages base)
  #:use-module (gnu packages cmake)
  #:use-module (gnu packages maths)
  #:use-module (gnu packages xml))

;(use-modules (guix git-download)
;             (guix build-system cmake)
;             (guix licenses)
;             (guix packages))

(define-public siconos
  (package
    (name "siconos")
    (version "4.3.x")
    (source (origin
              (method git-fetch)
              (uri (git-reference
                    (url "https://github.com/siconos/siconos")
                    (commit "HEAD")))
              (sha256 (base32
                       "0bgbz38lnxjhgghb5pfl0h1zjgly4vz5dc3wdr22550iym8clyab"))))
    (build-system cmake-build-system)
    (arguments
     '(#:build-type "Release"           ;Build without '-g' to save space.
                  #:configure-flags '()
                  #:tests? #f))                              ;XXX: no "test" target
    (inputs
     `(("swig" ,swig)
       ("boost" ,boost)
       ("gmp" ,gmp)
       ("openblas" ,openblas)
       ("lapack" ,lapack)
       ("python" ,python)
       ("python-numpy" ,python-numpy)
     ("python-scipy" ,python-scipy)
     ("gfortran" ,gfortran)
     ("gcc" ,gcc)
     ("gnu-make" ,gnu-make)
     ("cmake" ,cmake)))
    (home-page "https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos/index.html")
    (synopsis "Library for nonsmooth numerical simulation")
    (description
     "Siconos is an open-source scientific software primarily targeted at modeling and simulating nonsmooth dynamical systems in C++ and in Python:

    Mechanical systems (rigid or solid) with unilateral contact and Coulomb friction and impact (nonsmooth mechanics, contact dynamics, multibody systems dynamics or granular materials).
    Switched Electrical Circuit such as electrical circuits with ideal and piecewise linear components: power converter, rectifier, Phase-Locked Loop (PLL) or Analog-to-Digital converter.
    Sliding mode control systems.
    Biology (Gene regulatory network).

Other applications are found in Systems and Control (hybrid systems, differential inclusions, optimal control with state constraints), Optimization (Complementarity systems and Variational inequalities), Fluid Mechanics, and Computer Graphics.")
    (license "")))

