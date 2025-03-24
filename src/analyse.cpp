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

    cv::Mat kernel = cv::getStructuringElement(cv::MORPH_RECT, cv::Size(5, 5));
    cv::erode(mask, mask, kernel);


    cv::dilate(mask, mask, kernel);

    std::vector<std::vector<cv::Point>> contours;
    std::vector<cv::Vec4i> hierarchy;
    cv::findContours(mask, contours, hierarchy, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_SIMPLE);
    
    cv::Mat contourImage = cv::Mat(height, width, CV_8UC3, cv::Scalar(255, 255, 255));

    for (size_t i = 0; i < contours.size(); ++i) {

    cv::Mat pieceMask = cv::Mat::zeros(croppedImage.size(), CV_8UC1);
    cv::drawContours(pieceMask, contours, static_cast<int>(i), cv::Scalar(255), cv::FILLED);

    cv::Rect boundingBox = cv::boundingRect(contours[i]);

    cv::Mat piece;
    croppedImage.copyTo(piece, pieceMask); 

    cv::Mat croppedPiece = piece(boundingBox);
    cv::Mat croppedMask = pieceMask(boundingBox);

    croppedPiece.copyTo(contourImage(boundingBox), croppedMask);
}

int objIndex = 0;
for (const auto& contour : contours) {

    cv::Rect boundingBox = cv::boundingRect(contour);

    cv::Mat roi = contourImage(boundingBox);

    std::string filename = "data/object_" + std::to_string(objIndex) + ".png";
    cv::imwrite(filename, roi);
    objIndex++;
}

    cv::imshow("Image Originale (Recadrée)", croppedImage);
    cv::imshow("Masque basé sur la couleur dominante", mask);
    cv::imshow("Contours trouvés", contourImage);

    cv::waitKey(0);


    return 0;
}
