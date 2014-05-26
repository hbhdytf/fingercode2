% Filterbank-Based Fingerprint Matching (A.K.Jain, S.Prabhakar, L.Hong and S.Pankanti, 2000)
%
% Abstract
% With identity fraud in our society reaching unprecedented proportions and with an increasing emphasis on
% the emerging automatic personal identification applications, biometrics-based verification, especially
% fingerprint-based identification, is receiving a lot of attention. There are two major shortcomings
% of the traditional approaches to fingerprint representation. For a considerable fraction of population,
% the representations based on explicit detection of complete ridge structures in the fingerprint are
% difficult to extract automatically. The widely used minutiae-based representation does not utilize a
% significant component of the rich discriminatory information available in the fingerprints. Local ridge
% structures cannot be completely characterized by minutiae. Further, minutiae-based matching has difficulty
% in quickly matching two fingerprint images containing different number of unregistered minutiae points.
% The proposed filter-based algorithm uses a bank of Gabor filters to capture both local and global details
% in a fingerprint as a compact fixed length FingerCode. The fingerprint matching is based on the Euclidean
% distance between the two corresponding FingerCodes and hence is extremely fast. We are able to achieve a
% verification accuracy which is only marginally inferior to the best results of minutiae-based algorithms
% published in the open literature. Our system performs better than a state-of-the-art minutiae-based system
% when the performance requirement of the application system does not demand a very low false acceptance rate.
% Finally, we show that the matching performance can be improved by combining the decisions of the matchers
% based on complementary (minutiae-based and filter-based) fingerprint information.
%
% Index Terms: Biometrics, FingerCode, fingerprints, flow pattern, Gabor filters, matching, texture, verification.
%
% Type "fprec" on Matlab command window to start image processing.
% This source code provides a new, improved GUI respect to the previous release.
% Simulation parameters can be changed in this file:
%
% n_bands  : the number of concentric bands
% h_bands  : the width of each band in pixels
% n_arcs   : the number of arcs (each band has exactly n_arcs arcs)
% h_radius : the inner radius in pixel (the central band is not considered
%            because it is a too small area)
% num_disk : the number of Gabor filters
%
% n_sectors and h_lato rapresents respectively the total number of sectors
% (the length of the feature vector associated to each filter of filter-bank:
% there are num_disk Gabor filters) and the height of the cropped image in pixels.
% N_secors and h_lato should not be changed.
%
%
%   $Revision: 1.0 $  $Date: 2002.10.02  $
%   $Revision: 2.0 $  $Date: 2003.11.29  $
%   $Revision: 3.0 $  $Date: 2004.06.22  $
%   $Revision  4.0 $  $Date: 2005.03.21  $
%   $Revision  5.0 $  $Date: 2005.05.15  $
%   $Revision  5.1 $  $Date: 2005.05.21  $
%   $Revision  5.2 $  $Date: 2006.11.29  $   by Luigi Rosa
%                                            email:   luigi.rosa@tiscali.it
%                                            mobile:  +39 320 7214179
%                                            website: http://www.advancedsourcecode.com
%
%
%   Release           Date             Major features
%   5.2              2006.11.29       - Improved GUI
%                                     - Added PGM image format 
%                                     - Minor bug fixed
%   5.1              2005.05.21       - An improved alorithm used in fingerprint matching
%                                     - A generalized and optimized version of the algorithm for core point localization
%                                     - Improved database
%   5.0              2005.05.05       - The merging technique is greatly improved.
%                                       This new algorithm was tested on FVC2004 training fingerprint images
%                                     - Parameters for image segmentation are estimated automatically
%                                     - An improved algorithm is used when the rotated FingerCode
%                                       is added to database. Now this precedure does not introduce
%                                       any additional noise
%                                     - A faster fingerprint image acquisition when a new fingerprint image is added to database
%                                     - Implementation of recursive 1D and 2D recursive Gabor filtering
%                                     - Optimized pixel-wise orientation field estimation
%                                       ( 80% faster than Release 4.0 )
%                                     - List of fingerprint databases available on the web
%
%   4.0              2005.03.21       - An improved algorithm for core point detection,
%                                       based on a novel hybrid technique
%                                     - A better fingerprint segmentation
%                                       which makes use of morphological operations ((binary erosion and dilation)
%                                     - Exhaustive documentation of the implemented algorithms
%                                     - Improved GUI
%                                     - Better error management
%                                     - Image enhancement
%                                     - Orientation field estimation
%
%   3.0              2004.06.22       - Major bugs fixed
%                                     - New GUI
%                                     - Complex filtering techniques
%                                     - Improved core point determination
%                                     - Robustness against noise
%                                     - Modifiable simulation parameters
%
%   2.0              2003.11.29       - New GUI
%                                     - 8 Gabor filters 0 22.5 45 67.5 90 112.5 135 157.5 degrees
%                                     - Convolution is performed in frequency domain
%                                     - DataBase
%                                     - Fingerprint matching
%                                     - Error management
%
%
%
%
% Input fingerprint should be 256 x 256 image 8-bit grayscale @ 500 dpi.
% If these conditions are not verified some parameters in m-functions
% should be changed in a proper way (such as, for example, Gabor filter
% parameters in gabor2d_sub function). See the cited references for more
% details.
%
%   M-files included:
%
%     -fprec.m:                         This file, it initializes the entire image processing. The simulation
%                                       parameters can be changed in this main file.
%     -supercore7.m:                    A function which accepts an input image and determines the coordinates
%                                       of the core point. This implementation is based on a novel, improved algorithm.
%                                       The parameters for image segmentation are determined automatically.
%                                       This function is used when an input fingerprint image is added to database
%     -supercore7_list.m:               A function which accepts an input image and determines a list of the coordinates
%                                       of the core points candidates. This function is based on the algorithm
%                                       implemented in supercore7. This function is used when an input fingerprint image
%                                       is selected for fingerprint matching
%     -eliminacopie.m:                  This function, given an input sequence of N core points, returns M core points,
%                                       with M<= N , where all repetead core points are eliminated
%     -sogliavarianza.m                 Used to calculate the region of interest of input image
%     -filtraggiocomplesso.m            This function returns the complex filtering output
%     -mirror.m:                        A function which is used to "mirror" input image in  order to avoid undesired
%                                       boundary effects
%     -recrop.m:                        A function used to resize the mirrored filtered image
%     -conv2fft.m:                      This function performs 2D FFT-based convolution. Type "help conv2fft" on Matlab command
%                                       window for more details
%     -whichsector.m:                   A function used to determine (for each pixel of the cropped image) the corresponding
%                                       sectors of the concentric bands
%     -sector_norm.m:                   A function used to normalize input image and to calculate the features vector
%     -cropping.m:                      This function is used to cropp the input fingerprint image after the core point is
%                                       determinated
%     -gabor2d_sub.m:                   A function used to calculate the coefficients of the gabor 2D filters.
%     -vedicentro.m:                    This simple routines uses the M-function supercore7.m and it is used to display
%                                       the core point
%     -angular_filter_bank.m:           Precomputes angular filter bank
%     -compute_coherence.m:             Computes the coherence image
%     -fft_enhance_cubs.m:              Enhances the fingerprint image
%     -orientation_image_luigiopt.m:    Optimized pixel-wise orientation field estimation.
%     -orientation_image_luigiopt_var.m:Optimized pixel-wise orientation field estimation where the square moving window size can be chosen.
%                                       See orientation_image_luigiopt for more details.
%     -orientation_image_rao.m:         Block-wise orientation field estimation
%     -pseudo_matched_filter.m:         Implements root filtering technique to increase SNR of an image
%     -radial_filter_bank.m:            Precomputes angular filter bank
%     -smoothen_frequency_image.m:      Smoothens the frequency image through a process of diffusion
%     -smoothen_orientation_image.m:    Smoothens the orientation image through vectorial gaussian filtering
%     -view_orientation_image.m:        Displays the orientation image as a quiver plot
%     -gabor1d.m                        Recursive 1D Gabor filtering. Type "help gabor1d" on Matlab
%                                       command window for more details
%     -gabor2d.m                        Recursive 2D Gabor filtering. Type "help gabor2d" on Matlab
%                                       command window for more details
%     -gaborexamples.m                  Many examples of usage of recursive Gabor filtering
%
%
% A crucial step in fingerprint recognition is core point determination.
% If any error occurs while cropping image you can use the auxiliary m-file
% "vedicentro.m": it visualizes the input fingerprint and the core point
% calculated by the m-function "supercore7.m".
%
%
%  Notes:
%  The computational load can be significantly reduced by recursive filtering techniques.
%  For a complete publication list of Lucas J. van Vliet please visit the following URL:
%  http://www.ph.tn.tudelft.nl/~lucas/publications/papersLJvV.html
%  Here you will find articles concernings a recursive implementation of the Gaussian filter,
%  of the derivative Gaussian filter and of Gabor filter.
%
%  If you want to optimize the proposed method an excellent article is the following one:
%  Erian Bezhani, Dequn Sun, Jean-Luc Nagel, and Sergio Carrato, "Optimized filterbank
%  fingerprint recognition", Proc. SPIE Intern. Symp. Electronic Imaging 2003, 20-24
%  Jan. 2003, Santa Clara, California.
%
%  This code was developed using:
%    MATLAB Version 7.0.1   (R14) Service Pack 1
%    Operating System: Microsoft Windows 2000 Version 5.0 (Build 2195: Service Pack 4)
%    Java VM Version: Java 1.4.2_04 with Sun Microsystems Inc. Java HotSpot(TM) Client VM
%    Image Processing Toolbox  Version 5.0.1  (R14SP1) is required.
%
%  Please contribute if you find this software useful.
%  For the donation please visit my website:
%  http://utenti.lycos.it/matlab/fingerprint4.htm
%  Report bugs to luigi.rosa@tiscali.it
%
%
%   References:
%
%   Cheng Long Adam Wang, researcher
%   Fingerprint Recognition System
%   http://home.kimo.com.tw/carouse9/FRS.htm
%
%   A. K. Jain, S. Prabhakar, and S. Pankanti, "A Filterbank-based Representation for
%   Classification and Matching of Fingerprints", International Joint Conference on
%   Neural Networks (IJCNN), pp. 3284-3285, Washington DC, July 10-16, 1999.
%   http://www.cse.msu.edu/~prabhaka/publications.html
%
%   "Fingerprint Classification and Matching Using a Filterbank", Salil Prabhakar
%   A DISSERTATION Submitted to Michigan State University in partial fulfillment
%   of the requirements for the degree of DOCTOR OF PHILOSOPHY, Computer
%   Science & Engineering, 2001
%   http://biometrics.cse.msu.edu/SalilThesis.pdf
%
%   Final Report 18-551 (Spring 1999) Fingerprint Recognition Group Number 19
%   Markus Adhiwiyogo, Samuel Chong, Joseph Huang, Weechoon Teo
%   http://www.ece.cmu.edu/~ee551/Old_projects/projects/s99_19/finalreport.html
%
%   Kenneth Nilsson and Josef Bigun, "Localization of corresponding points in
%   fingerprints by complex filtering", Pattern Recognition Letters, 24 (2003) 2135-2144
%   School of Information Science, Computer and Electrical Engineering (IDE), Halmstad
%   University, P.O. Box 823, SE-301 18, Halmstad, Sweden.
%   http://www.hh.se/staff/josef/public/publications/nilsson03prl.pdf
%
%   S. Chikkerur,C. Wu and Govindaraju, "A Systematic approach for feature extraction
%   in fingerprint images", ICBA 2004
%   http://www.eng.buffalo.edu/~ssc5/research.html#fingerprints
%   See also Chikkerur's submission at Matlab File Exchange
%   http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=6830&objectType=file
%
%   I.T. Young, L.J. van Vliet, M. van Ginkel, Recursive Gabor filtering, in: A. Sanfeliu, J.J. Villanueva,
%   M. Vanrell, R. Alquezar, T. Huang, J. Serra (eds.), ICPR15, Proc. 15th Int. Conference on Pattern Recognition
%   (Barcelona, Spain, Sep.3-7), vol. 3, Image, Speech, and Signal Processing, IEEE Computer Society Press,
%   Los Alamitos, 2000, 342-345.
%   http://www.ph.tn.tudelft.nl/~lucas/publications/2000/ICPR2000TYLVMG/ICPR2000TYLVMG.html
%
%   I.T. Young, L.J. van Vliet, M. van Ginkel, Recursive Gabor filtering, IEEE Transactions on Signal
%   Processing, Vol. 50, No. 11, 2798-2805, 2002.
%   http://www.ph.tn.tudelft.nl/~lucas/publications/2002/IEEETSP2002TYLVMG/IEEETSP2002TYLVMG.html
%
% ************************************************************************
% *  This is my postal address:                                          *
% *  Luigi Rosa                                                          *
% *  Via Centrale 35                                                     *
% *  67042 Civita di Bagno                                               *
% *  L'Aquila --- ITALY                                                  *
% *                                                                      *
% *  mobile +39 320 7214179                                              *
% *  email luigi.rosa@tiscali.it                                         *
% *  website http://www.advancedsourcecode.com                           *
% ************************************************************************
%
%
%
