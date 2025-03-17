/***************************************************************************************************
 * @file  analyse.cpp
 * @brief Contains the main program of the project
 **************************************************************************************************/
 #include <opencv4/opencv2/opencv.hpp>
 #include <iostream>
 
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

int main() {

    cv::Mat inputImage = cv::imread("data/puzzle.jpg");
    if (inputImage.empty()) {
        std::cerr << "Erreur lors du chargement de l'image." << std::endl;
        return -1;
    }

    cv::Mat imageHSV;
    cv::cvtColor(inputImage, imageHSV, cv::COLOR_BGR2HSV);

    //cv::Mat blurred;
    //cv::GaussianBlur(imageHSV, blurred, cv::Size(0, 0), 0);  

    cv::Scalar dominantColor = getDominantColor(imageHSV);
    std::cout << "Teinte dominante : " << dominantColor[0] << std::endl;


    int hueRange = 5; 
    cv::Scalar lowerBound(dominantColor[0] - 2 , 170, 170);
    cv::Scalar upperBound(dominantColor[0] + 5, 255, 255);

    cv::Mat mask;
    cv::inRange(imageHSV, lowerBound, upperBound, mask);

    // Affichage
    cv::imshow("Image Originale", imageHSV);
    cv::imshow("Masque basé sur la couleur dominante", mask);
    cv::waitKey(0);

    return 0;
}

