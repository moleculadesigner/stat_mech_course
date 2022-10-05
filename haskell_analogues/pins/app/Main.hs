module Main where

import           System.Random (RandomGen, mkStdGen, uniformR, randomRs)
import           Data.List (unfoldr, sort)
import           Options.Applicative
import           Data.Time.Clock



main :: IO ()
main = do
    time <- getCurrentTime >>= pure . floor . utctDayTime
    gen <- pure $ mkStdGen time
    options <- execParser opts
    let left_ = left options
    let right_ = right options
    let sigma_ = sigma options
    let nPins_ = nPins options
    let nRuns_ = nRuns options
    print 
        $ take nRuns_
        $ validStates left_ right_ sigma_
        $ states gen nPins_ left_ right_


pack :: Int -> [a] -> [[a]]
pack n = unfoldr (\l -> Just (take n l, drop n l))


states :: RandomGen g => g -> Int -> Double -> Double -> [[Double]]
states gen n left right = pack n $ randomRs (left, right) gen


checkState :: Double -> Double -> Double -> [Double] -> Bool
checkState start end sigma space = and $ differences start (sort space)
  where
    differences st [] = []
    differences st [pin] = [((pin - st) >= 2*sigma) && ((end - pin) >= 2*sigma)]
    differences st (pin:pins) =
        ((pin - st) >= 2*sigma)
        : (differences pin pins)


validStates :: Double -> Double -> Double -> [[Double]] -> [[Double]]
validStates left right sigma = 
    filter $ checkState (left - sigma) (right + sigma) sigma


-- Argument parser
data Parameters = Parameters
    { left  :: Double
    , right :: Double
    , sigma :: Double
    , nPins :: Int
    , nRuns :: Int
    }

leftP :: Parser Double
leftP = option auto
    ( long "left"
    <> short 'l'
    <> metavar "L"
    <> help "Left bound of the rope"
    )

rightP :: Parser Double
rightP = option auto
    ( long "right"
    <> short 'r'
    <> metavar "R"
    <> help "Right bound of the rope"
    )

sigmaP :: Parser Double
sigmaP = option auto
    ( long "sigma"
    <> short 's'
    <> metavar "S"
    <> help "Pin width"
    )

nPinsP :: Parser Int
nPinsP = option auto
    ( long "pins"
    <> short 'p'
    <> metavar "P"
    <> help "Pin number"
    )

nRunsP :: Parser Int
nRunsP = option auto
    ( long "runs"
    <> short 'n'
    <> metavar "N"
    <> help "Number of successful placements"
    )

parameters :: Parser Parameters
parameters = Parameters
    <$> leftP
    <*> rightP
    <*> sigmaP
    <*> nPinsP
    <*> nRunsP

opts :: ParserInfo Parameters
opts = info (parameters <**> helper)
  ( fullDesc
  <> progDesc "Sample pins on a rope"
  <> header "Specify simulation parameters" )