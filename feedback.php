<?php
// feedback.php
// PHP version of feedback.cgi (anonymous feedback form)

error_reporting(E_ALL);
ini_set("display_errors", 1);

$name      = "Lars Winther Christensen";
$fromName  = "Feedback Form";
$email     = "lars.w.christensen@ttu.edu";

$subject   = "Math 2360 student feedback";

if ($_SERVER["REQUEST_METHOD"] === "POST" && isset($_POST['button']) && $_POST['button'] === "Submit") {
    $feedback = trim($_POST['feedback']);
    if (!empty($feedback)) {
        $message = nl2br(htmlspecialchars($feedback));
        $headers  = "From: $fromName <mathdept@ttu.edu>\r\n";
        $headers .= "Reply-To: mathdept@ttu.edu\r\n";
        $headers .= "MIME-Version: 1.0\r\n";
        $headers .= "Content-Type: text/html; charset=UTF-8\r\n";

        if (mail($email, $subject, $message, $headers)) {
            echo "<!DOCTYPE html><html><head><title>Email Sent</title></head><body>
                  <h2>Thanks for the feedback</h2>
                  <p>Your message was successfully sent to <i>$name</i>.</p>
                  </body></html>";
        } else {
            echo "<p>⚠️ Error: Email could not be sent. Please contact the administrator.</p>";
        }
    } else {
        echo "<p>Please enter feedback before submitting.</p>";
    }
} else {
    echo "<!DOCTYPE html><html><head><title>Anonymous Feedback</title></head><body>
          <h2>Welcome to the anonymous feedback page</h2>
          <p>Write your comments and suggestions in the form below.</p>
          <p>Press <i>Submit</i> to send an anonymous email to $name.</p>
          <form method='POST' action=''>
              <textarea name='feedback' rows='8' cols='72'></textarea><br><br>
              <input type='submit' name='button' value='Submit'>
          </form></body></html>";
}
?>

