function send_Email()
%--------------------------------------------------------------------------
% Benachrichtigungs E-Mail
%--------------------------------------------------------------------------
disp('==============================================================');
disp('                        Sending E-Mail                        ');
disp('==============================================================');
% Define these variables appropriately:
mail = 'gmail.com'; %Your GMail email address
password = ''; %Your GMail password
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','E_mail',mail);
setpref('Internet','SMTP_Username',mail);
setpref('Internet','SMTP_Password',password);
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');
% Send the email. Note that the first input is the address you are sending the email to
sendmail('tommy-c22452@web.de',(['MATLAB Job: 123 für ' num2str(123) '^3 und delta = ' num2str(123)]),(['splashingalgorithm completed successfully: grid size ', num2str(123), ' run in ', num2str(123)]))
%--------------------------------------------------------------------------
end 
