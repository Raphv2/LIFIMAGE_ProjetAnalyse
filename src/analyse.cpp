/***************************************************************************************************
 * @file  analyse.cpp
 * @brief Contains the main program of the project
 **************************************************************************************************/
 #include <opencv4/opencv2/opencv.hpp>
 #include <iostream>
 #include <vector>
#include <cstdlib> 
#include <ctime>   
 
 cv::Scalar getDominantColor(cv::Mat& imageHSV) {
     std::vector<cv::Mat> hsvChannels;
     cv::split(imageHSV, hsvChannels);
 
     int histSize = 180; 
     float range[] = {0, 180}; 
     const float* histRange = {range};
     cv::Mat hist;
 
     cv::calcHist(&hsvChannels[0], 1, 0, cv::Mat(), hist, 1, &histSize, &histRange, true, false);
 
     double maxVal = 0;
     cv::Point maxLoc; 
     cv::minMaxLoc(hist, 0, &maxVal, 0, &maxLoc);
     
     int dominantHue = maxLoc.y; 
 
     return cv::Scalar(dominantHue, 100, 100); 
 }

 cv::Scalar getRandomColor() {
    return cv::Scalar(std::rand() % 256, std::rand() % 256, std::rand() % 256);
}
 
 int main() {
     cv::Mat inputImage = cv::imread("data/puzzle.jpg");
     if (inputImage.empty()) {
         std::cerr << "Erreur lors du chargement de l'image." << std::endl;
         return -1;
     }

     int width = inputImage.cols;
     int height = inputImage.rows;
 
     cv::Rect resize(10 , 10 , width * 0.9, height - 20); 
 
     cv::Mat croppedImage = inputImage(resize);
 
     cv::Mat imageHSV;
     cv::cvtColor(croppedImage, imageHSV, cv::COLOR_BGR2HSV);
 
     cv::Scalar dominantColor = getDominantColor(imageHSV);
     std::cout << "Teinte dominante : " << dominantColor[0] << std::endl;
 
     int hueRange = 5; 
     cv::Scalar lowerBound(dominantColor[0] - hueRange, 170, 170);
     cv::Scalar upperBound(dominantColor[0] + hueRange, 255, 255);
 
     cv::Mat mask;
     cv::inRange(imageHSV, lowerBound, upperBound, mask);
     cv::bitwise_not(mask, mask);

     std::vector<std::vector<cv::Point>> contours;
     std::vector<cv::Vec4i> hierarchy;
     cv::findContours(mask, contours, hierarchy, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_SIMPLE);
 
     cv::Mat contourImage = cv::Mat(height, width, CV_8UC3, cv::Scalar(255, 255, 255));
    
     for (size_t i = 0; i < contours.size(); ++i) {
        cv::Scalar color = getRandomColor(); 
        cv::drawContours(contourImage, contours, static_cast<int>(i), color, cv::FILLED); 
    }

    cv::Mat kernel = cv::getStructuringElement(cv::MORPH_RECT, cv::Size(5, 5));
    cv::erode(contourImage, contourImage, kernel);


    cv::dilate(contourImage, contourImage, kernel);

    
     // Affichage
     cv::imshow("Image Originale (Recadrée)", croppedImage);
     cv::imshow("Masque basé sur la couleur dominante", mask);
     cv::imshow("Contours trouvés", contourImage);

 
     cv::waitKey(0);

 
     return 0;
 }
 